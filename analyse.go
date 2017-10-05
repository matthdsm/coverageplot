package coverageplot

import (
	"encoding/json"
	"fmt"
	"github.com/biogo/hts/sam"
	"github.com/brentp/irelate/parsers"
	"github.com/brentp/xopen"
	log "github.com/sirupsen/logrus"
	"io"
	"math"
	"strconv"
	"strings"
	"time"
)

type ip struct {
	chrom string
	start uint32
	end   uint32
}

func (p ip) Chrom() string {
	return p.chrom
}
func (p ip) Start() uint32 {
	return p.start
}
func (p ip) End() uint32 {
	return p.end
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

type readinterval struct {
	chrom string
	start int
	end   int
	first bool
}

func Analyse(bamFile string, intervalsFile string, format string) {
	//get intervals
	log.Info(fmt.Sprintf("Opening file intervals: %q", intervalsFile))
	rdr, err := xopen.Ropen(intervalsFile)
	check(err)
	defer rdr.Close()
	var intervals []Interval
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
		var infos []Info
		switch format {
		case "bed":
			infos = nil
		case "gff":
			err = json.Unmarshal([]byte(toks[3]), &infos)
			if err != nil {
				fmt.Println("error:", err)
			}
		default:
			panic(fmt.Sprintf("Format %s is not valid.", format))
		}
		//log.Debug( fmt.Sprintf("%+v\n", infos))
		intervals = append(intervals, Interval{toks[0], s, e, infos})

	}

	// open bam file
	b, _ := parsers.NewBamQueryable(bamFile) // 0-based half open
	intervalCount := len(intervals)
	intervalProgress := -1
	template := "[%s%s] %d%% %d/%d %s"
	startTime := time.Now()
	log.Info(fmt.Sprintf("Going to process %d intervals", intervalCount))
	for intervalPos, i := range intervals {
		log.Debug(fmt.Sprintf("Calculating coverage for interval %d: %s", intervalPos, i.String()))
		reg := ip{i.chrom, uint32(i.start), uint32(i.end)}
		processReads(b, reg)
		intervalPercentage := int((float64(intervalPos) / float64(intervalCount)) * 100)
		if intervalPercentage >= intervalProgress+1 {
			progressSize := int(math.Ceil(float64(intervalPercentage) * 0.75))
			progressBar := strings.Repeat("#", progressSize)
			restSize := 75 - progressSize
			restBar := strings.Repeat("-", restSize)

			log.Info(fmt.Sprintf(template, progressBar, restBar, intervalPercentage, intervalPos, intervalCount, time.Since(startTime)))
			intervalProgress = intervalProgress + 1
		}
	}

}

func processReads(bamFile *parsers.BamQueryable, interval ip) {
	intervalSize := interval.end - interval.start
	var CurrCoverage []int
	CurrCoverage = make([]int, intervalSize)

	q, err := bamFile.Query(interval)
	if err != nil {
		log.Warn("Region did not contain any reads")
	}
	heap := make(map[string]readinterval)

	j := 0
	for rec, e := q.Next(); e == nil; rec, e = q.Next() {
		read := rec.(*parsers.Bam)

		if read.Flags&sam.Unmapped != 0 {
			log.Debug(fmt.Sprintf("skipped, read is unmappped %d", read.Flags&sam.Duplicate))
			continue
		}

		if read.Flags&sam.Secondary != 0 {
			log.Debug(fmt.Sprintf("skipped, read is Secondary %d", read.Flags&sam.Duplicate))
			continue
		}

		if read.Record.MapQ < 1 {
			log.Debug("skipped, read mapping is low quality")
			continue
		}

		if read.Flags&sam.Duplicate != 0 {
			log.Debug(fmt.Sprintf("skipped, read is duplicate %d", read.Flags&sam.Duplicate))
			continue
		}

		// Skip if not a proper pair
		if read.Flags&sam.Paired != 0 &&
			read.Flags&sam.ProperPair != 0 &&
			read.Flags&sam.MateUnmapped != 0 {
			log.Debug(fmt.Sprintf("skipped, read badPair check, flag %v", read.Flags))
			continue
		}

		// Skip if wrong coordinates
		if ((read.Pos < read.MatePos) && read.Flags&sam.Reverse != 0) ||
			((read.MatePos < read.Pos) && read.Flags&sam.MateReverse != 0) {
			log.Debug(fmt.Sprintf("skipped, chemically designed: left on positive strand, right on reverse one, flag %v", read.Flags))
			//chemically designed: left on positive strand, right on reverse one
			continue
		}

		//log.Debug( fmt.Sprintf("Found read on %s:%d-%d with MQ %d and flags %v\n", rec.Chrom(), rec.Start(), rec.End(), read.MapQ, read.Flags,))
		end := int(read.End())

		overlap := (2 * read.Len()) - abs(read.TempLen)
		//duplicate := 1
		//if read.Flags&sam.Duplicate == 0 {
		//	duplicate = 0
		//}
		//if read.MapQ <= 40 {
		//	log.Debug( fmt.Sprintf("Found the other read with lower MQ: %s\n", read.MapQ))
		//}
		log.Debug(fmt.Sprintf("%s\t%v\t%s\t%d\t%s:%d-%d\t%d\tmate\t%s:%d\t%d\t%d\t%d",
			read.Name,
			read.Flags,
			read.Cigar,
			read.Record.MapQ,
			read.Ref.Name(),
			read.Pos,
			end,
			end-read.Pos,
			read.MateRef.Name(),
			read.MatePos,
			read.TempLen,
			overlap,
			end-read.MatePos,
		//read.AuxFields,
		))
		if read.Flags&sam.Read1 != 0 {
			if readi, ok := heap[read.Name]; ok {
				if !readi.first {
					log.Debug(fmt.Sprintf("Found the other read2 on the heap, checking overlap for read: %s", read.Name))
					//add coverage minus overlap
					CurrCoverage = addcoverage(CurrCoverage, removeoverlap(readi, readinterval{read.Ref.Name(), read.Pos, end, true}), interval)
					delete(heap, read.Name)
				}
			} else {
				heap[read.Name] = readinterval{read.Ref.Name(), read.Pos, end, true}
				log.Debug(fmt.Sprintf("Added first read %s to heap", read.Name))
			}
			CurrCoverage = addcoverage(CurrCoverage, readinterval{read.Ref.Name(), read.Pos, end, true}, interval)

		}
		if read.Flags&sam.Read2 != 0 {
			if readi, ok := heap[read.Name]; ok {
				log.Debug(fmt.Sprintf("Found the other read1 on the heap, checking overlap for read: %s", read.Name))
				//add coverage minus overlap
				CurrCoverage = addcoverage(CurrCoverage, removeoverlap(readinterval{read.Ref.Name(), read.Pos, end, false}, readi), interval)
				delete(heap, read.Name)
			} else {
				heap[read.Name] = readinterval{read.Ref.Name(), read.Pos, end, false}
				log.Debug(fmt.Sprintf("Added second read %s to heap", read.Name))
			}
		}
		j++
	}
	log.Debug(fmt.Sprintf("%d reads left on the heap", len(heap)))
	// add coverage from reads on the heap
	for _, heapread := range heap {
		CurrCoverage = addcoverage(CurrCoverage, heapread, interval)
	}
	log.Debug(fmt.Sprintf("%d reads found in the interval", j))
	log.Debug(fmt.Sprintf("Coverage string: %v", CurrCoverage))
	for i, coverage_i := range CurrCoverage {
		fmt.Printf("%s\t%d\t%d\tNA\t%d\n", interval.chrom, i+int(interval.start), 1+i+int(interval.start), coverage_i)
	}
}

func removeoverlap(reada readinterval, readb readinterval) readinterval {
	// remove overlap from reada/read2
	if reada.first || !readb.first {
		log.Debug("Can't substract read 2 from read 1")
		return reada
	}
	if reada.chrom != readb.chrom {
		log.Debug("Reads in pair on different chromosomes")
		return reada
	}
	overlap := 0
	if reada.start >= readb.start { // read b is leftmost
		// remove overlap from the start of a
		overlap = readb.end - reada.start
		if overlap > 0 {
			reada.start = readb.end
		}
	} else { // read a is leftmost
		// remove overlap from the end of a
		overlap = reada.end - readb.start
		if overlap > 0 {
			reada.end = readb.start
		}

	}
	log.Debug(fmt.Sprintf("%d overlap", overlap))
	//readinterval{read.Ref.Name(), read.Pos, end, true}
	//	chrom string
	// start int
	// end   int
	// first bool
	//if reada.start <
	return reada
}

func addcoverage(CurrCoverage []int, read readinterval, interval ip) []int {
	//readinterval{read.Ref.Name(), read.Pos, end, true}
	log.Debug(fmt.Sprintf("Add coverage %s:%d..%d", read.chrom, read.start, read.end))
	coverageNetStart := max(read.start, int(interval.start)) - int(interval.start)
	coverageNetEnd := min(read.end, int(interval.end)) - int(interval.start)
	log.Debug(fmt.Sprintf("net coverage %d..%d", coverageNetStart+int(interval.start), coverageNetEnd+int(interval.start)))

	for i := coverageNetStart; i < coverageNetEnd; i++ {
		CurrCoverage[i]++
	}
	return CurrCoverage
}

//func resetChromCoverage() {
//    _currChromName = "";
//    _currChromSize = 0 ;
//    std::vector<DEPTH>().swap(_currChromCoverage);
//}
// void BedGenomeCoverage::AddBlockedCoverage(const vector<BED> &bedBlocks) {
//     vector<BED>::const_iterator bedItr = bedBlocks.begin();
//     vector<BED>::const_iterator bedEnd = bedBlocks.end();
//     for (; bedItr != bedEnd; ++bedItr) {
//         // the end - 1 must be done because BamAncillary::getBamBlocks
//         // returns ends uncorrected for the genomeCoverageBed data structure.
//         // ugly, but necessary.
//         AddCoverage(bedItr->start, bedItr->end - 1);
//     }
// }

// void BedGenomeCoverage::AddCoverage(int start, int end) {
//     // process the first line for this chromosome.
//     // make sure the coordinates fit within the chrom
//     if (start < _currChromSize)
//         _currChromCoverage[start].starts++;
//     if (end >= 0 && end < _currChromSize)
//         _currChromCoverage[end].ends++;
//     else
//         _currChromCoverage[_currChromSize-1].ends++;
// }

//const (
//	Paired        uint = iota // The read is paired in sequencing, no matter whether it is mapped in a pair.
//	ProperPair                // The read is mapped in a proper pair.
//	Unmapped                  // The read itself is unmapped; conflictive with ProperPair.
//	MateUnmapped              // The mate is unmapped.
//	Reverse                   // The read is mapped to the reverse strand.
//	MateReverse               // The mate is mapped to the reverse strand.
//	Read1                     // This is read1.
//	Read2                     // This is read2.
//	Secondary                 // Not primary alignment.
//	QCFail                    // QC failure.
//	Duplicate                 // Optical or PCR duplicate.
//	Supplementary             // Supplementary alignment, indicates alignment is part of a chimeric alignment.
//)
