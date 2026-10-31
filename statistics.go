package coverageplot

import (
	"bufio"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"log/slog"
	"os/exec"
	"strconv"
	"strings"
)

var defaultBucketLabels = []float64{0, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}

var execCommand = exec.Command

// BedInfo contains exon and transcript metadata.
type BedInfo struct {
	Gene         string
	RefSeqId     string
	TranscriptId string
	ExonId       string
	Rank         string
}

// BedInterval represents a genomic range in 0-based half-open coordinates.
type BedInterval struct {
	chrom string
	start int
	end   int
	info  []BedInfo
}

// Size returns the length in base pairs of the interval.
func (i BedInterval) Size() int {
	return i.end - i.start
}

// Statistic contains calculated coverage metrics for a region or total panel.
type Statistic struct {
	Gene_panel           string
	Genome_build         string
	Chromosome           string
	Start                int
	End                  int
	Gene                 string
	Transcript_ensembl   string
	Transcript_refseq    string
	Exon_ensemble        string
	Rank                 string
	Size                 int
	Coverage             []float64
	Min                  float64
	Max                  float64
	Mean                 float64
	StdDev               float64
	Quartiles            Quartiles
	Zero_bases           int
	CummulativeHistogram []int
	Histogram            []int
}

// StatisticOutput is the top-level container for total and regional statistics.
type StatisticOutput struct {
	Name    string
	Total   Statistic
	Regions []Statistic
}

func summarizeStatistic(coverage []float64, histo *Histogram) Statistic {
	var stat Statistic
	stat.Size = len(coverage)
	stat.Quartiles, _ = Quartile(coverage)
	stat.Mean, _ = Mean(coverage)
	stat.Min, _ = Min(coverage)
	stat.Max, _ = Max(coverage)
	stat.StdDev, _ = StandardDeviation(coverage)
	stat.CummulativeHistogram = histo.BucketsCumm
	stat.Histogram = histo.Buckets
	if len(histo.Buckets) > 0 {
		stat.Zero_bases = histo.Buckets[0]
	}
	return stat
}

// Statistics executes coverage calculation and outputs JSON to stdout.
func Statistics(intervalsFile string, coverageBed string, panelName string, total bool) error {
	var panelStats StatisticOutput
	if total && intervalsFile == "" {
		res, err := CalculateTotalStatistics(coverageBed, panelName)
		if err != nil {
			return err
		}
		panelStats = res
	} else {
		intervals, err := ParseRegions(intervalsFile)
		if err != nil {
			return err
		}
		res, err := CalculateStatistics(intervals, intervalsFile, coverageBed, panelName, total)
		if err != nil {
			return err
		}
		panelStats = res
	}

	data, err := json.Marshal(panelStats)
	if err != nil {
		return fmt.Errorf("marshal statistics JSON: %w", err)
	}
	fmt.Println(string(data))
	return nil
}

// ParseRegions reads and parses intervals from a panel BED file.
func ParseRegions(intervalsFile string) ([]BedInterval, error) {
	slog.Debug("Opening intervals bed file", "file", intervalsFile)
	rdr, err := openReader(intervalsFile)
	if err != nil {
		return nil, fmt.Errorf("open intervals file %q: %w", intervalsFile, err)
	}
	defer rdr.Close()

	var intervals []BedInterval
	br := bufio.NewReader(rdr)
	for {
		line, readErr := br.ReadString('\n')
		trimmed := strings.TrimRight(line, "\r\n")
		if trimmed != "" {
			toks := strings.Split(trimmed, "\t")
			if len(toks) >= 3 {
				s, err1 := strconv.Atoi(toks[1])
				e, err2 := strconv.Atoi(toks[2])
				if err1 == nil && err2 == nil {
					if e < s {
						slog.Warn("Region start must be leftmost coordinate; swapping", "start", s, "end", e)
						s, e = e, s
					}

					var infos []BedInfo
					if len(toks) > 3 {
						for _, infoValues := range strings.Split(toks[3], "|") {
							infoArr := strings.Split(infoValues, ";")
							gene := infoArr[0]
							if len(infoArr) > 4 {
								infos = append(infos, BedInfo{
									Gene:         gene,
									RefSeqId:     infoArr[1],
									TranscriptId: infoArr[2],
									ExonId:       infoArr[3],
									Rank:         infoArr[4],
								})
							} else if len(infoArr) > 0 {
								infos = append(infos, BedInfo{Gene: infoArr[0]})
							}
						}
					}
					c := strings.TrimPrefix(toks[0], "chr")
					intervals = append(intervals, BedInterval{chrom: c, start: s, end: e, info: infos})
				}
			}
		}

		if readErr != nil {
			if errors.Is(readErr, io.EOF) {
				break
			}
			return nil, fmt.Errorf("read intervals file: %w", readErr)
		}
	}
	return intervals, nil
}

// CalculateStatistics calculates regional and summary statistics by querying tabix on coverageBed.
func CalculateStatistics(intervals []BedInterval, intervalsFile string, coverageBed string, panelName string, total bool) (StatisticOutput, error) {
	slog.Debug("Calculating statistics", "coverage", coverageBed, "intervals", intervalsFile)

	var panelStats StatisticOutput
	panelStats.Name = panelName
	panelStats.Total.Gene_panel = panelName
	panelStats.Total.Genome_build = "GRCh38"

	if len(intervals) == 0 {
		return panelStats, nil
	}

	var cmd *exec.Cmd
	if intervalsFile == "intervals" {
		region := fmt.Sprintf("%s:%d-%d", intervals[0].chrom, intervals[0].start, intervals[0].end)
		slog.Debug("Running tabix for region", "region", region)
		cmd = execCommand("tabix", coverageBed, region)
	} else {
		slog.Debug("Running tabix for intervals file", "file", intervalsFile)
		cmd = execCommand("tabix", "-R", intervalsFile, coverageBed)
	}

	stdout, err := cmd.StdoutPipe()
	if err != nil {
		return panelStats, fmt.Errorf("tabix stdout pipe: %w", err)
	}
	if err := cmd.Start(); err != nil {
		return panelStats, fmt.Errorf("start tabix: %w", err)
	}
	defer func() {
		_ = cmd.Wait()
	}()

	rdrCov := bufio.NewReader(stdout)
	var totalCoverage []float64
	var lastInterval BedInterval
	totalHisto := NewHistogram(defaultBucketLabels)

	for _, interval := range intervals {
		currCoverage := make([]float64, interval.Size())
		isFirstSegment := true

		for {
			covLine, readErr := rdrCov.ReadString('\n')
			trimmed := strings.TrimRight(covLine, "\r\n")
			if readErr != nil && errors.Is(readErr, io.EOF) && trimmed == "" {
				break
			}

			covToks := strings.Split(trimmed, "\t")
			if len(covToks) < 4 {
				if readErr != nil {
					break
				}
				continue
			}

			start, err1 := strconv.Atoi(covToks[1])
			end, err2 := strconv.Atoi(covToks[2])
			if err1 != nil || err2 != nil {
				continue
			}

			if isFirstSegment && (interval.start < start || interval.start >= end) {
				continue
			}
			isFirstSegment = false

			coverage, errParse := strconv.ParseFloat(covToks[3], 64)
			if errParse != nil {
				continue
			}

			for i := start; i < end; i++ {
				idx := i - interval.start
				if idx >= 0 && idx < len(currCoverage) {
					currCoverage[idx] = coverage
				}
			}

			if end >= interval.end {
				break
			}
			if readErr != nil {
				break
			}
		}

		if !total {
			histo := NewHistogram(defaultBucketLabels)
			histo.AddAll(currCoverage)
			stat := summarizeStatistic(currCoverage, histo)
			stat.Gene_panel = panelName
			stat.Genome_build = "GRCh38"
			stat.Chromosome = interval.chrom
			stat.Coverage = currCoverage
			stat.Start = interval.start
			stat.End = interval.end
			if len(interval.info) > 0 {
				stat.Gene = interval.info[0].Gene
				stat.Transcript_ensembl = interval.info[0].TranscriptId
				stat.Transcript_refseq = interval.info[0].RefSeqId
				stat.Exon_ensemble = interval.info[0].ExonId
				stat.Rank = interval.info[0].Rank
			}
			panelStats.Regions = append(panelStats.Regions, stat)
		}

		if interval.chrom == lastInterval.chrom && interval.start < lastInterval.end {
			indexCorrected := lastInterval.end - interval.start
			if indexCorrected > interval.Size() {
				continue
			}
			if indexCorrected > 0 {
				totalCoverage = append(totalCoverage, currCoverage[indexCorrected:]...)
				lastInterval = interval
				totalHisto.AddAll(currCoverage[indexCorrected:])
			}
		} else {
			totalCoverage = append(totalCoverage, currCoverage...)
			lastInterval = interval
			totalHisto.AddAll(currCoverage)
		}
	}

	panelStats.Total = summarizeStatistic(totalCoverage, totalHisto)
	panelStats.Total.Gene_panel = panelName
	panelStats.Total.Genome_build = "GRCh38"
	panelStats.Total.Start = 0
	panelStats.Total.End = len(totalCoverage)

	return panelStats, nil
}

// CalculateTotalStatistics computes whole-file coverage statistics in a streaming single pass.
func CalculateTotalStatistics(coverageBed string, panelName string) (StatisticOutput, error) {
	slog.Debug("Calculating total statistics stream", "coverage", coverageBed)
	rdr, err := openReader(coverageBed)
	if err != nil {
		return StatisticOutput{}, fmt.Errorf("open coverage bed %q: %w", coverageBed, err)
	}
	defer rdr.Close()

	var runningStats Stats
	runningStats.Reset()
	q := NewTargeted(0.25, 0.50, 0.75)
	totalHisto := NewHistogram(defaultBucketLabels)

	br := bufio.NewReader(rdr)
	for {
		covLine, readErr := br.ReadString('\n')
		trimmed := strings.TrimRight(covLine, "\r\n")
		if trimmed != "" {
			covToks := strings.Split(trimmed, "\t")
			if len(covToks) >= 4 {
				start, err1 := strconv.Atoi(covToks[1])
				end, err2 := strconv.Atoi(covToks[2])
				count, err3 := strconv.ParseFloat(covToks[3], 64)
				if err1 == nil && err2 == nil && err3 == nil {
					span := end - start
					if span > 0 {
						totalHisto.AddCount(count, span)
						runningStats.AddCount(count, uint64(span))
						q.AddCount(count, uint64(span))
					}
				}
			}
		}

		if readErr != nil {
			if errors.Is(readErr, io.EOF) {
				break
			}
			return StatisticOutput{}, fmt.Errorf("read coverage bed: %w", readErr)
		}
	}

	var returnStats Statistic
	returnStats.Gene_panel = panelName
	returnStats.Genome_build = "GRCh38"
	returnStats.Start = 0
	returnStats.End = int(runningStats.Count())
	returnStats.Size = int(runningStats.Count())
	returnStats.Mean = runningStats.Mean()
	returnStats.Min = runningStats.Min()
	returnStats.Max = runningStats.Max()
	returnStats.StdDev = runningStats.Stddev()
	returnStats.Quartiles = Quartiles{
		Q1: q.Query(0.25),
		Q2: q.Query(0.50),
		Q3: q.Query(0.75),
	}
	returnStats.CummulativeHistogram = totalHisto.BucketsCumm
	returnStats.Histogram = totalHisto.Buckets
	if len(totalHisto.Buckets) > 0 {
		returnStats.Zero_bases = totalHisto.Buckets[0]
	}

	return StatisticOutput{
		Name:  panelName,
		Total: returnStats,
	}, nil
}
