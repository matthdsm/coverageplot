package coverageplot

import (
	"bufio"
	"fmt"
	"github.com/brentp/xopen"
	"github.com/montanaflynn/stats"
	log "github.com/sirupsen/logrus"
	"io"
	"os/exec"
	"strconv"
	"strings"
)

type BedInfo struct { // NCF1;NM_000265.5;ENST00000289473;ENSE00003510291;8
	Gene         string
	RefSeqId     string
	TranscriptId string
	ExonId       string
	Rank         string
}

type BedInterval struct {
	chrom string
	start int
	end   int
	info  []BedInfo
}

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
	Quartiles            stats.Quartiles
	Zero_bases           int
	CummulativeHistogram []int
	Histogram            []int
}

func Statistics(intervalsFile string, coverageBed string, panelName string) {
	//get intervals
	log.Info(fmt.Sprintf("Opening file intervals bed: %q", intervalsFile))
	rdr, err := xopen.Ropen(intervalsFile)
	check(err)
	defer rdr.Close()

	log.Info(fmt.Sprintf("Opening file coverage bed: %q", coverageBed))
	log.Debug(fmt.Sprintf("Using command: |tabix -R %s %s", intervalsFile, coverageBed))
	//fmt.Sprintf("|tabix %s %s", intervalsFile, coverageBed)
	//rdr_cov, err := xopen.Ropen("|head -n50 data/D1706395.chr_5000lines.coverage.bed")
	var rdr_test io.Reader
	var cmd *exec.Cmd
	cmd = exec.Command("sh", "-c", fmt.Sprintf("tabix -R %s %s", intervalsFile, coverageBed))
	rdr_test, err = cmd.StdoutPipe()
	err = cmd.Start()
	check(err)
	rdr_cov := bufio.NewReader(rdr_test)
	//check(err)
	//defer rdr_cov.Close()

	//var intervals []BedInterval
	var totLength int = 0
	var TotalCoverage []float64

	for {
		line, err := rdr.ReadString('\n')
		if err == io.EOF {
			break
		} else {
			check(err)
		}
		toks := strings.Split(strings.TrimSuffix(line, "\n"), "\t")
		s, err := strconv.Atoi(toks[1])
		e, err := strconv.Atoi(toks[2])
		if e < s {
			log.Warn("region start has to be the left most coordinate, but a switched region found, fix your input file")
			t := s
			s = e
			e = t
		}
		var infos []BedInfo
		var gene string
		for i, info_values := range strings.Split(toks[3], "|") {
			info_values_arr := strings.Split(info_values, ";")
			if i == 0 {
				gene = info_values_arr[0]
			}
			infos = append(infos, BedInfo{gene, info_values_arr[1], info_values_arr[2], info_values_arr[3], info_values_arr[4]})
		}
		c := strings.TrimPrefix(toks[0], "chr")
		interval := BedInterval{c, s, e, infos}
		totLength = totLength + interval.Size()
		var CurrCoverage []float64
		CurrCoverage = make([]float64, interval.Size())
		log.Debug(fmt.Sprintf("Create coverage for %dbp region", interval.Size()))
		for {
			cov_line, err := rdr_cov.ReadString(byte('\n'))
			cov_toks := strings.Split(strings.TrimSuffix(cov_line, "\n"), "\t")
			if err == io.EOF {
				log.Debug("Break EOF")
				break
			} else {
				check(err)
			}
			start, err := strconv.Atoi(cov_toks[1])
			rel_start := start - interval.start
			count, err := strconv.ParseFloat(cov_toks[4], 64)
			CurrCoverage[rel_start] = count
			if start == (interval.end - 1) {
				break
			}
		}
		var stat Statistic
		stat.Gene_panel = panelName
		stat.Genome_build = "GRCh38"
		stat.Chromosome = interval.chrom
		stat.Coverage = CurrCoverage
		stat.Start = interval.start
		stat.End = interval.end
		stat.Gene = interval.info[0].Gene
		stat.Transcript_ensembl = interval.info[0].TranscriptId
		stat.Transcript_refseq = interval.info[0].RefSeqId
		stat.Exon_ensemble = interval.info[0].ExonId
		stat.Rank = interval.info[0].Rank
		stat.Size = interval.Size()
		stat.Quartiles, _ = stats.Quartile(CurrCoverage)
		stat.Mean, _ = stats.Mean(CurrCoverage)
		stat.Min, _ = stats.Min(CurrCoverage)
		stat.Max, _ = stats.Max(CurrCoverage)
		stat.StdDev, _ = stats.StandardDeviation(CurrCoverage)
		var histo Histogram
		histo.init([]float64{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100})
		histo.AddAll(CurrCoverage)
		stat.CummulativeHistogram = histo.BucketsCumm
		stat.Histogram = histo.Buckets
		stat.Zero_bases = stat.Histogram[0]
		if stat.Zero_bases > 0 {
			log.Warn(fmt.Sprintf("Zero coverage for %d bases, exon %s, %s:%d-%d", stat.Zero_bases, stat.Exon_ensemble, stat.Chromosome, stat.Start, stat.End))
		}

		// my @my_header = ("Gene_panel","Genome_build","Chromosome", "start", "end", "Gene", "Transcripts_ensembl", "Transcripts_refseq", "Exon_ensemble", "Rank", "Target_region_length", "Min_coverage", "Max_coverage", "Mean_coverage", "Median_coverage", "Zero_coverage_bases", "Proportion_covered","Target_region_length_with_coverage_above_20", "Percentage_with_coverage_above_20", "Target_region_length_with_coverage_above_30","Percentage_with_coverage_above_30"); # to output
		log.Debug(fmt.Sprintf("Stats: %v", stat))
		TotalCoverage = append(TotalCoverage, CurrCoverage...)

	}
	// total of all intervals
	log.Info(fmt.Sprintf("Total selected size (in kbp): %.3f", float64(totLength)/1000))
	var totalStat Statistic
	totalStat.Gene_panel = panelName
	totalStat.Genome_build = "GRCh38"
	totalStat.Start = 0
	totalStat.End = len(TotalCoverage)
	totalStat.Size = totLength
	totalStat.Quartiles, _ = stats.Quartile(TotalCoverage)
	totalStat.Mean, _ = stats.Mean(TotalCoverage)
	totalStat.Min, _ = stats.Min(TotalCoverage)
	totalStat.Max, _ = stats.Max(TotalCoverage)
	totalStat.StdDev, _ = stats.StandardDeviation(TotalCoverage)
	var totalHisto Histogram
	totalHisto.init([]float64{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100})
	totalHisto.AddAll(TotalCoverage)
	totalStat.CummulativeHistogram = totalHisto.BucketsCumm
	totalStat.Histogram = totalHisto.Buckets
	totalStat.Zero_bases = totalStat.Histogram[0]
	log.Debug(fmt.Sprintf("Total Stats: %v", totalStat))

}

func (i BedInterval) Size() int {
	return (i.end - i.start)
}
