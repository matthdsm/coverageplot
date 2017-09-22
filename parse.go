package coverageplot

import (
	"encoding/json"
	"fmt"
	"github.com/brentp/xopen"
	log "github.com/sirupsen/logrus"
	"io"
	"regexp"
	"strconv"
	"strings"
)

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func Parse(gffFile string) {
	log.Info(fmt.Sprintf("Opening file: %q", gffFile))
	for intervals := range exons(gffFile) {
		for _, i := range intervals {
			fmt.Printf("%s\n", i.String())
		}
	}
}

func Merge(intervalFile string, expand int, format string) {
	log.Info(fmt.Sprintf("Opening file: %q", intervalFile))
	log.Info(fmt.Sprintf("Writing output to stdout in %s format", format))
	log.Info(fmt.Sprintf("Expanding exons with %d bp at both ends", expand))
	rdr, err := xopen.Ropen(intervalFile)
	check(err)
	defer rdr.Close()
	var totLength int = 0
	var chrom string
	var start int = 0
	var end int = 4294967295
	var info []string
	var strand []string
	for {
		line, err := rdr.ReadString('\n')
		if err == io.EOF {
			break
		} else {
			check(err)
		}
		toks := strings.Split(strings.TrimSuffix(line, "\n"), "\t")
		if len(toks) < 4 {
			toks = append(toks, "")
		}
		s, err := strconv.Atoi(toks[1])
		e, err := strconv.Atoi(toks[2])
		if e < s {
			log.Warn("region start has to be the left most coordinate, but a switched region found, fix your input file")
		}
		if end == 4294967295 {
			chrom = strings.TrimPrefix(toks[0], "chr")
			start = s
			end = e
			info = nil
			info = append(info, toks[3])
			strand = nil
			if len(toks) > 5 {
				strand = append(strand, toks[5])
			} else {
				strand = append(strand, "+")
			}
		}
		if chrom == strings.TrimPrefix(toks[0], "chr") && s >= start && s <= (end+expand+expand) {
			start = min(s, start)
			if end == 4294967295 {
				end = e
			} else {
				end = max(end, e)
			}
			if info != nil && toks[3] != info[len(info)-1] {
				info = append(info, toks[3])
			}
			if len(toks) > 5 {
				strand = append(strand, toks[5])
			}
		} else {
			fmt.Print(FormatLine(format, chrom, start-expand, end+expand, info, "0", strand))
			totLength = totLength + (end - start)
			chrom = strings.TrimPrefix(toks[0], "chr")
			start = s
			end = e
			info = nil
			info = append(info, toks[3])
			strand = nil
			if len(toks) > 5 {
				strand = append(strand, toks[5])
			} else {
				strand = append(strand, "+")
			}
		}
	}
	//print last interval
	fmt.Print(FormatLine(format, chrom, start-expand, end+expand, info, "0", strand))
	totLength = totLength + (end - start)
	// total of all intervals
	log.Info(fmt.Sprintf("Total size (in kbp): %d", totLength/1000))

}

func FormatLine(format string, chrom string, start int, end int, info []string, score string, strand []string) string {
	switch format {
	case "bed":
		// if bed the coordinates are already 0-based
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%s\t%s\n", chrom, max(0, start), end, strings.Join(info, "|"), score, strings.Join(strand, ","))
	case "gff":
		// if gff input convert from 1-based coordinates to 0-based coordinates
		return fmt.Sprintf("%s\t%d\t%d\t[%s]\t%s\t%s\n", chrom, max(0, start-1), end, strings.Join(info, ","), score, strand[0])
	default:
		panic(fmt.Sprintf("Format %s is not valid.\n", format))
	}
}

type Info struct {
	Gene          string
	Product       string
	Transcript_id string
	Gbkey         string
	Id            string // Dbxref=GeneID:100288667,HGNC:HGNC:37756
	Hgnc          string // Dbxref=GeneID:100288667,HGNC:HGNC:37756
	Genbank       string // Dbxref=GeneID:100133331,Genbank:NR_028327.1
}

type Interval struct {
	chrom string
	start int
	end   int
	info  []Info
}

func (i Interval) String() string {
	info, _ := json.Marshal(i.info)
	return fmt.Sprintf("%s\t%d\t%d\t%s", i.chrom, i.start, i.end, string(info))
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func exons(gffFile string) chan []Interval {
	ch := make(chan []Interval)
	nIntervals := 50

	go func() {
		rdr, err := xopen.Ropen(gffFile)
		check(err)
		defer rdr.Close()
		regions := make([]Interval, 0, nIntervals)
		for {
			line, err := rdr.ReadString('\n')
			if err == io.EOF {
				break
			} else {
				check(err)
			}
			toks := strings.Split(strings.TrimSuffix(line, "\n"), "\t")
			if toks[0][0] == '#' {
				continue
			}
			re_chrom := regexp.MustCompile(`NC_0+([1-9][0-9]*)\.\d+`)
			parsed_chrom := re_chrom.FindStringSubmatch(toks[0])
			if len(parsed_chrom) < 1 {
				continue
			}
			if parsed_chrom[1] == "12920" { //chrMT
				continue
			}
			c := parsed_chrom[1]

			s, err := strconv.Atoi(toks[3])
			e, err := strconv.Atoi(toks[4])

			var info Info
			// parse info - col 8
			for _, info_tok := range strings.Split(toks[8], ";") {
				info_kv := strings.Split(info_tok, "=")
				switch info_kv[0] {
				case "gene":
					info.Gene = info_kv[1]
				case "product":
					info.Product = info_kv[1]
				case "transcript_id":
					info.Transcript_id = info_kv[1]
				case "gbkey":
					info.Gbkey = info_kv[1]
				case "Dbxref":
					for _, xref_tok := range strings.Split(info_kv[1], ",") {
						xref_kv := strings.Split(xref_tok, ":")
						if xref_kv[0] == "GeneID" {
							info.Id = xref_kv[1]
						}
						if xref_kv[0] == "HGNC" {
							info.Hgnc = xref_kv[len(xref_kv)-1]
						}
						if xref_kv[0] == "Genbank" {
							info.Genbank = xref_kv[1]
						}

					}
				default:
					continue
				}
			}
			infos := []Info{info}
			var r Interval
			if s > e {
				r = Interval{c, e, s, infos}
			} else {
				r = Interval{c, s, e, infos}
			}

			regions = append(regions, r)
			if len(regions) == nIntervals {
				ch <- regions
				regions = make([]Interval, 0, nIntervals)
			}

		}
		if len(regions) > 0 {
			ch <- regions
		}
		close(ch)

	}()
	return ch

}
