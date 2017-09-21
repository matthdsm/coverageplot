package coverageplot

import (
	"encoding/json"
	"fmt"
	"github.com/biogo/hts/sam"
	"github.com/brentp/irelate/parsers"
	"github.com/brentp/xopen"
	"io"
	"os"
	"strconv"
	"strings"
)

var CurrCoverage []int
var CurrRegionStart int
var CurrRegionEnd int

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

func Analyse(bamFile string, intervalsFile string) {
	//get intervals
	fmt.Fprintf(os.Stderr, "Opening file intervals: %q\n", intervalsFile)
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
		err = json.Unmarshal([]byte(toks[3]), &infos)
		if err != nil {
			fmt.Println("error:", err)
		}
		//fmt.Fprintf(os.Stderr, "%+v\n", infos)
		intervals = append(intervals, Interval{toks[0], s, e, infos})

	}

	// open bam file
	b, _ := parsers.NewBamQueryable(bamFile) // 0-based half open
	for _, i := range intervals {
		fmt.Printf("### Calculating coverage for:\n%s\n", i.String())
		reg := ip{i.chrom, uint32(i.start), uint32(i.end)}
		processReads(b, reg)
	}

}

func processReads(bamFile *parsers.BamQueryable, interval ip) {
	intervalSize := interval.end - interval.start
	CurrCoverage = make([]int, intervalSize)
	q, err := bamFile.Query(interval)
	if err != nil {
		fmt.Println(os.Stderr, "Region did not contain any reads\n")
	}
	heap := make(map[string]readinterval)

	j := 0
	for rec, e := q.Next(); e == nil; rec, e = q.Next() {
		read := rec.(*parsers.Bam)

		if read.Flags&sam.Unmapped != 0 {
			fmt.Fprintf(os.Stderr, "skipped, read is unmappped %d\n", read.Flags&sam.Duplicate)
			continue
		}

		if read.Flags&sam.Secondary != 0 {
			fmt.Fprintf(os.Stderr, "skipped, read is Secondary %d\n", read.Flags&sam.Duplicate)
			continue
		}

		if read.Record.MapQ < 1 {
			fmt.Fprint(os.Stderr, "skipped, read mapping is low quality\n")
			continue
		}

		//if read.Flags&sam.Duplicate != 0 {
		//	fmt.Fprintf(os.Stderr, "skipped, read is duplicate %d\n", read.Flags&sam.Duplicate)
		//	continue
		//}

		// Skip if not a proper pair
		if read.Flags&sam.Paired != 0 &&
			read.Flags&sam.ProperPair != 0 &&
			read.Flags&sam.MateUnmapped != 0 {
			fmt.Fprintf(os.Stderr, "skipped, read badPair check, flag %v\n", read.Flags)
			continue
		}

		// Skip if wrong coordinates
		if ((read.Pos < read.MatePos) && read.Flags&sam.Reverse != 0) ||
			((read.MatePos < read.Pos) && read.Flags&sam.MateReverse != 0) {
			fmt.Fprintf(os.Stderr, "skipped, chemically designed: left on positive strand, right on reverse one, flag %v\n", read.Flags)
			//chemically designed: left on positive strand, right on reverse one
			continue
		}

		//fmt.Fprintf(os.Stderr, "Found read on %s:%d-%d with MQ %d and flags %v\n", rec.Chrom(), rec.Start(), rec.End(), read.MapQ, read.Flags,)
		end := int(read.End())
		overlap := (2 * read.Len()) - abs(read.TempLen)
		//duplicate := 1
		//if read.Flags&sam.Duplicate == 0 {
		//	duplicate = 0
		//}
		//if read.MapQ <= 40 {
		//	fmt.Fprintf(os.Stderr, "Found the other read with lower MQ: %s\n", read.MapQ)
		//}
		fmt.Printf("%s\t%v\t%s\t%d\t%s:%d-%d\t%d\tmate\t%s:%d\t%d\t%d\t%d\n",
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
		)
		if read.Flags&sam.Read1 != 0 {
			if readi, ok := heap[read.Name]; ok {
				if !readi.first {
					fmt.Fprintf(os.Stderr, "Found the other read2 on the heap, checking overlap for read: %s\n", read.Name)
					//add coverage minus overlap
					addcoverage(removeoverlap(readi, readinterval{read.Ref.Name(), read.Pos, end, true}), interval)
					delete(heap, read.Name)
				}
			} else {
				heap[read.Name] = readinterval{read.Ref.Name(), read.Pos, end, true}
				fmt.Fprintf(os.Stderr, "Added first read %s to heap\n", read.Name)
			}
			addcoverage(readinterval{read.Ref.Name(), read.Pos, end, true}, interval)

		}
		if read.Flags&sam.Read2 != 0 {
			if readi, ok := heap[read.Name]; ok {
				fmt.Fprintf(os.Stderr, "Found the other read1 on the heap, checking overlap for read: %s\n", read.Name)
				//add coverage minus overlap
				addcoverage(removeoverlap(readinterval{read.Ref.Name(), read.Pos, end, false}, readi), interval)
				delete(heap, read.Name)
			} else {
				heap[read.Name] = readinterval{read.Ref.Name(), read.Pos, end, false}
				fmt.Fprintf(os.Stderr, "Added second read %s to heap\n", read.Name)
			}
		}
		j++
	}
	fmt.Fprintf(os.Stderr, "%d reads left on the heap\n", len(heap))
	// add coverage from reads on the heap
	for _, heapread := range heap {
		addcoverage(heapread, interval)
	}
	fmt.Fprintf(os.Stderr, "%d reads found in the interval\n", j)
	fmt.Fprintf(os.Stderr, "Coverage string: %v\n", CurrCoverage)
	//for i, coverage_i := range CurrCoverage {
	//	fmt.Fprintf(os.Stderr, "chr1\t%d\t%d\n", i+CurrRegionStart, coverage_i)
	//}

}

func removeoverlap(reada readinterval, readb readinterval) readinterval {
	// remove overlap from reada/read2
	if reada.first || !readb.first {
		fmt.Fprint(os.Stderr, "Can't substract read 2 from read 1\n")
		return reada
	}
	if reada.chrom != readb.chrom {
		fmt.Fprint(os.Stderr, "Reads in pair on different chromosomes\n")
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
	fmt.Fprintf(os.Stderr, "%d overlap\n", overlap)
	//readinterval{read.Ref.Name(), read.Pos, end, true}
	//	chrom string
	// start int
	// end   int
	// first bool
	//if reada.start <
	return reada
}

func addcoverage(read readinterval, interval ip) {
	//readinterval{read.Ref.Name(), read.Pos, end, true}
	fmt.Fprintf(os.Stderr, "Add coverage %s:%d..%d\n", read.chrom, read.start, read.end)
	coverageNetStart := max(read.start, int(interval.start)) - int(interval.start)
	coverageNetEnd := min(read.end, int(interval.end)) - int(interval.start)
	fmt.Fprintf(os.Stderr, "net coverage %d..%d\n", coverageNetStart+int(interval.start), coverageNetEnd+int(interval.start))

	for i := coverageNetStart; i < coverageNetEnd; i++ {
		CurrCoverage[i]++
	}
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
