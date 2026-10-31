package coverageplot

import (
	"flag"
	"log/slog"
	"math"
	"os"
	"os/exec"
	"slices"
	"strings"
	"testing"
)

var debugTest bool

func TestMain(m *testing.M) {
	flag.BoolVar(&debugTest, "debug", false, "Enable debug logging for tests")
	flag.Parse()

	logLevel := slog.LevelInfo
	if debugTest {
		logLevel = slog.LevelDebug
	}
	slog.SetDefault(slog.New(slog.NewTextHandler(os.Stderr, &slog.HandlerOptions{Level: logLevel})))

	os.Exit(m.Run())
}

func mockExecCommand(command string, args ...string) *exec.Cmd {
	cs := []string{"-test.run=TestHelperProcess", "--", command}
	cs = append(cs, args...)
	cmd := exec.Command(os.Args[0], cs...)

	testType := "standard"
	for _, arg := range args {
		if strings.Contains(arg, "identical_end_test.bed.gz") {
			testType = "identical_end"
			break
		} else if strings.Contains(arg, "overlap_test.bed.gz") {
			testType = "overlap"
			break
		} else if strings.Contains(arg, "loop_termination_test.bed.gz") || strings.Contains(arg, "loop_termination.bed.gz") {
			testType = "loop_termination"
			break
		}
	}
	cmd.Env = []string{"GO_WANT_HELPER_PROCESS=1", "TEST_TYPE=" + testType}
	return cmd
}

func TestHelperProcess(t *testing.T) {
	if os.Getenv("GO_WANT_HELPER_PROCESS") != "1" {
		return
	}

	testType := os.Getenv("TEST_TYPE")

	switch testType {
	case "identical_end":
		mockOutput := `chr1	90	100	15
chr1	100	110	35
chr1	110	130	30
chr1	130	150	55
chr1	110	130	30
chr1	130	150	55
`
		os.Stdout.WriteString(mockOutput)
	case "overlap":
		mockOutput := `chr1	90	100	10
chr1	100	110	30
chr1	110	120	40
chr1	120	150	50
chr1	120	150	50
`
		os.Stdout.WriteString(mockOutput)
	case "real_data":
		data, err := os.ReadFile("data/TabixOutput_exome.bed")
		if err == nil {
			os.Stdout.Write(data)
		}
	case "loop_termination":
		mockOutput := `chr1	90	95	5
chr1	95	100	10
chr1	100	110	20
chr1	110	120	30
chr1	120	125	40
chr2	95	100	10
chr2	100	120	40
chr2	120	120	70
chr2	120	125	60
chr2	125	125	70
`
		os.Stdout.WriteString(mockOutput)
	case "tango2_data":
		data, err := os.ReadFile("data/tabixOut_TANGO2_test.exome.bed")
		if err == nil {
			os.Stdout.Write(data)
		} else {
			mockOutput := `chr22	42016645	42016670	68
chr22	42016670	42016700	72
chr22	42016700	42016730	65
chr22	42016730	42016770	80
chr22	42018460	42018490	90
chr22	42018490	42018520	95
chr22	42018520	42018560	105
chr22	42019817	42019840	115
chr22	42019840	42019870	120
`
			os.Stdout.WriteString(mockOutput)
		}
	default:
		mockOutput := `chr1	90	100	10
chr1	100	110	30
chr1	110	120	40
chr1	120	130	50
chr1	200	210	20
chr1	210	220	25
chr1	220	230	30
chr1	230	232	30
chr2	100	110	45
chr2	113	120	50
chr2	120	130	55
chr2	130	140	50
`
		os.Stdout.WriteString(mockOutput)
	}
	os.Exit(0)
}

func assertEqual[T comparable](t *testing.T, expected, actual T, msg ...string) {
	t.Helper()
	if expected != actual {
		t.Fatalf("expected %v, got %v: %s", expected, actual, strings.Join(msg, " "))
	}
}

func assertSliceEqual[T comparable](t *testing.T, expected, actual []T) {
	t.Helper()
	if !slices.Equal(expected, actual) {
		t.Fatalf("expected slice %v, got %v", expected, actual)
	}
}

func assertInEpsilon(t *testing.T, expected, actual, epsilon float64) {
	t.Helper()
	diff := math.Abs(expected - actual)
	if diff > epsilon {
		t.Fatalf("expected %v, got %v (diff: %v > eps: %v)", expected, actual, diff, epsilon)
	}
}

func assertTrue(t *testing.T, cond bool, msg string) {
	t.Helper()
	if !cond {
		t.Fatalf("expected true: %s", msg)
	}
}

func TestCalculateStatistics(t *testing.T) {
	origExecCommand := execCommand
	execCommand = mockExecCommand
	defer func() { execCommand = origExecCommand }()

	tmpfile, err := os.CreateTemp("", "test.bed")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmpfile.Name())

	intervals := []BedInterval{
		{
			chrom: "chr1",
			start: 95,
			end:   125,
			info:  []BedInfo{{Gene: "GENE1", TranscriptId: "ENST01", RefSeqId: "NM_001", ExonId: "ENSE01", Rank: "1"}},
		},
		{
			chrom: "chr1",
			start: 201,
			end:   231,
			info:  []BedInfo{{Gene: "GENE2", TranscriptId: "ENST02", RefSeqId: "NM_002", ExonId: "ENSE02", Rank: "1"}},
		},
		{
			chrom: "chr2",
			start: 100,
			end:   130,
			info:  []BedInfo{{Gene: "GENE3", TranscriptId: "ENST03", RefSeqId: "NM_003", ExonId: "ENSE03", Rank: "1"}},
		},
	}

	result, err := CalculateStatistics(intervals, tmpfile.Name(), "test.bed.gz", "TestPanel", false)
	if err != nil {
		t.Fatalf("CalculateStatistics failed: %v", err)
	}

	assertSliceEqual(t, []float64{10, 10, 10, 10, 10, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 50, 50, 50, 50, 50}, result.Regions[0].Coverage)
	assertEqual(t, "TestPanel", result.Name)
	assertEqual(t, 3, len(result.Regions))
	assertEqual(t, "GENE1", result.Regions[0].Gene)
	assertEqual(t, 95, result.Regions[0].Start)
	assertEqual(t, 125, result.Regions[0].End)
	assertEqual(t, 30, result.Regions[0].Size)
	assertEqual(t, 10.0, result.Regions[0].Min)
	assertEqual(t, 50.0, result.Regions[0].Max)
	assertInEpsilon(t, 33.33, result.Regions[0].Mean, 0.01)
	assertSliceEqual(t, []int{0, 0, 5, 0, 10, 10, 5, 0, 0, 0, 0, 0}, result.Regions[0].Histogram)
	assertSliceEqual(t, []int{0, 0, 5, 5, 15, 25, 30, 30, 30, 30, 30, 30}, result.Regions[0].CummulativeHistogram)

	assertSliceEqual(t, []float64{20, 20, 20, 20, 20, 20, 20, 20, 20, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30}, result.Regions[1].Coverage)
	assertEqual(t, "GENE2", result.Regions[1].Gene)
	assertEqual(t, 201, result.Regions[1].Start)
	assertEqual(t, 231, result.Regions[1].End)
	assertEqual(t, 30, result.Regions[1].Size)
	assertEqual(t, 20.0, result.Regions[1].Min)
	assertEqual(t, 30.0, result.Regions[1].Max)
	assertInEpsilon(t, 25.33, result.Regions[1].Mean, 0.01)

	assertSliceEqual(t, []float64{45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 0, 0, 0, 50, 50, 50, 50, 50, 50, 50, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55}, result.Regions[2].Coverage)
	assertEqual(t, "GENE3", result.Regions[2].Gene)
	assertEqual(t, 100, result.Regions[2].Start)
	assertEqual(t, 130, result.Regions[2].End)
	assertEqual(t, 30, result.Regions[2].Size)
	assertEqual(t, 0.0, result.Regions[2].Min)
	assertEqual(t, 55.0, result.Regions[2].Max)
	assertInEpsilon(t, 45.0, result.Regions[2].Mean, 0.01)

	assertEqual(t, 90, result.Total.Size)
	assertEqual(t, 0.0, result.Total.Min)
	assertEqual(t, 55.0, result.Total.Max)

	totalResult, err := CalculateStatistics(intervals, tmpfile.Name(), "test.bed.gz", "TestPanel", true)
	if err != nil {
		t.Fatalf("CalculateStatistics total failed: %v", err)
	}
	assertEqual(t, 0, len(totalResult.Regions))
	assertEqual(t, 90, totalResult.Total.Size)
}

func TestCalculateStatisticsWithOverlappingIntervals(t *testing.T) {
	origExecCommand := execCommand
	execCommand = mockExecCommand
	defer func() { execCommand = origExecCommand }()

	tmpfile, err := os.CreateTemp("", "test.bed")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmpfile.Name())

	intervals := []BedInterval{
		{
			chrom: "chr1",
			start: 95,
			end:   125,
			info:  []BedInfo{{Gene: "GENE1", TranscriptId: "ENST01", RefSeqId: "NM_001", ExonId: "ENSE01", Rank: "1"}},
		},
		{
			chrom: "chr1",
			start: 120,
			end:   150,
			info:  []BedInfo{{Gene: "GENE2", TranscriptId: "ENST02", RefSeqId: "NM_002", ExonId: "ENSE02", Rank: "1"}},
		},
	}

	result, err := CalculateStatistics(intervals, tmpfile.Name(), "overlap_test.bed.gz", "TestPanel", false)
	if err != nil {
		t.Fatalf("CalculateStatistics failed: %v", err)
	}

	assertEqual(t, 2, len(result.Regions))
	assertEqual(t, 55, result.Total.Size)
	assertSliceEqual(t, []float64{10, 10, 10, 10, 10, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 50, 50, 50, 50, 50}, result.Regions[0].Coverage)
	assertEqual(t, 50.0, result.Regions[1].Min)
	assertEqual(t, 50.0, result.Regions[1].Max)
}

func TestCalculateStatisticsWithOverlappingIdenticalEndIntervals(t *testing.T) {
	origExecCommand := execCommand
	execCommand = mockExecCommand
	defer func() { execCommand = origExecCommand }()

	tmpfile, err := os.CreateTemp("", "test.bed")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmpfile.Name())

	intervals := []BedInterval{
		{
			chrom: "chr1",
			start: 95,
			end:   150,
			info:  []BedInfo{{Gene: "GENE1", TranscriptId: "ENST01", RefSeqId: "NM_001", ExonId: "ENSE01", Rank: "1"}},
		},
		{
			chrom: "chr1",
			start: 120,
			end:   150,
			info:  []BedInfo{{Gene: "GENE2", TranscriptId: "ENST02", RefSeqId: "NM_002", ExonId: "ENSE02", Rank: "1"}},
		},
	}

	result, err := CalculateStatistics(intervals, tmpfile.Name(), "identical_end_test.bed.gz", "TestPanel", false)
	if err != nil {
		t.Fatalf("CalculateStatistics failed: %v", err)
	}

	assertSliceEqual(t, []float64{30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55}, result.Regions[1].Coverage)
	assertEqual(t, 2, len(result.Regions))
	assertEqual(t, 55, result.Total.Size)
	assertEqual(t, 55.0, result.Total.Max)
	assertEqual(t, 15.0, result.Total.Min)
}

func TestCalculateStatisticsWithEdgeCases(t *testing.T) {
	emptyResult, err := CalculateStatistics([]BedInterval{}, "intervals", "test.bed.gz", "EmptyPanel", false)
	if err != nil {
		t.Fatalf("unexpected error for empty intervals: %v", err)
	}
	assertEqual(t, 0, len(emptyResult.Regions))
	assertEqual(t, 0, emptyResult.Total.Size)
}

func TestCalculateStatisticsWithCriticalLoopTermination(t *testing.T) {
	origExecCommand := execCommand
	defer func() { execCommand = origExecCommand }()

	execCommand = func(command string, args ...string) *exec.Cmd {
		cs := []string{"-test.run=TestHelperProcess", "--", command}
		cs = append(cs, args...)
		cmd := exec.Command(os.Args[0], cs...)
		cmd.Env = []string{"GO_WANT_HELPER_PROCESS=1", "TEST_TYPE=loop_termination"}
		return cmd
	}

	tmpfile, err := os.CreateTemp("", "test.bed")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmpfile.Name())

	intervals := []BedInterval{
		{
			chrom: "chr1",
			start: 95,
			end:   125,
			info:  []BedInfo{{Gene: "BOUNDARY1", TranscriptId: "ENST01", RefSeqId: "NM_001", ExonId: "ENSE01", Rank: "1"}},
		},
		{
			chrom: "chr2",
			start: 95,
			end:   125,
			info:  []BedInfo{{Gene: "BOUNDARY2", TranscriptId: "ENST02", RefSeqId: "NM_002", ExonId: "ENSE02", Rank: "1"}},
		},
	}

	result, err := CalculateStatistics(intervals, tmpfile.Name(), "loop_termination.bed.gz", "TerminationTest", false)
	if err != nil {
		t.Fatalf("CalculateStatistics failed: %v", err)
	}

	assertEqual(t, 10.0, result.Regions[0].Coverage[0])
	assertEqual(t, 40.0, result.Regions[0].Coverage[29])
	assertEqual(t, 40.0, result.Regions[0].Max)
	assertEqual(t, 60.0, result.Regions[1].Coverage[29])

	assertEqual(t, 2, len(result.Regions))
	assertEqual(t, 60, result.Total.Size)
	assertEqual(t, 60.0, result.Total.Max)
	assertEqual(t, 10.0, result.Total.Min)
}

func TestCalculateStatisticsWithRealData(t *testing.T) {
	intervalsFile := "data/intervals_per_exon.bed"
	if _, err := os.Stat(intervalsFile); os.IsNotExist(err) {
		t.Skip("Skipping test: intervals file not found")
	}

	origExecCommand := execCommand
	defer func() { execCommand = origExecCommand }()

	execCommand = func(command string, args ...string) *exec.Cmd {
		cs := []string{"-test.run=TestHelperProcess", "--", command}
		cs = append(cs, args...)
		cmd := exec.Command(os.Args[0], cs...)
		cmd.Env = []string{"GO_WANT_HELPER_PROCESS=1", "TEST_TYPE=real_data"}
		return cmd
	}

	intervals, err := ParseRegions(intervalsFile)
	if err != nil {
		t.Fatalf("ParseRegions failed: %v", err)
	}
	if len(intervals) == 0 {
		t.Fatal("No intervals parsed from the real intervals file")
	}

	result, err := CalculateStatistics(intervals, intervalsFile, "real_data_test.bed.gz", "RealDataTest", false)
	if err != nil {
		t.Fatalf("CalculateStatistics failed: %v", err)
	}

	assertEqual(t, "RealDataTest", result.Name)
	assertTrue(t, len(result.Regions) > 0, "Should have at least one region")
	assertTrue(t, result.Total.Size > 0, "Total size should be greater than 0")
	assertEqual(t, 5.0, result.Total.Min)
	assertEqual(t, 259.0, result.Total.Max)
	assertInEpsilon(t, 97.72, result.Total.Mean, 0.01)

	for _, region := range result.Regions[:5] {
		assertTrue(t, region.End > region.Start, "End should be greater than start")
		assertEqual(t, region.End-region.Start, region.Size, "Size should match end-start")
		assertEqual(t, len(region.Coverage), region.Size, "Coverage array length should match size")
	}
}

func TestCalculateStatisticsWithTango2Data(t *testing.T) {
	intervalsFile := "data/TANGO2_intervals_per_exon.bed"
	if _, err := os.Stat(intervalsFile); os.IsNotExist(err) {
		t.Skip("Skipping test: TANGO2 intervals file not found")
	}

	origExecCommand := execCommand
	defer func() { execCommand = origExecCommand }()

	execCommand = func(command string, args ...string) *exec.Cmd {
		cs := []string{"-test.run=TestHelperProcess", "--", command}
		cs = append(cs, args...)
		cmd := exec.Command(os.Args[0], cs...)
		cmd.Env = []string{"GO_WANT_HELPER_PROCESS=1", "TEST_TYPE=tango2_data"}
		return cmd
	}

	intervals, err := ParseRegions(intervalsFile)
	if err != nil {
		t.Fatalf("ParseRegions failed: %v", err)
	}
	if len(intervals) == 0 {
		t.Fatal("No intervals parsed from the TANGO2 intervals file")
	}

	result, err := CalculateStatistics(intervals, intervalsFile, "data/tabixOut_TANGO2_test.exome.bed", "TANGO2Test", false)
	if err != nil {
		t.Fatalf("CalculateStatistics failed: %v", err)
	}

	assertEqual(t, "TANGO2Test", result.Name)
	assertTrue(t, len(result.Regions) > 0, "Should have at least one region")
	assertTrue(t, result.Total.Size > 0, "Total size should be greater than 0")

	assertInEpsilon(t, 73.8, result.Regions[3].Mean, 0.01)
	assertInEpsilon(t, 77.23, result.Regions[4].Mean, 0.01)
	assertInEpsilon(t, 108.41, result.Regions[5].Mean, 0.01)

	for _, region := range result.Regions[3:6] {
		assertTrue(t, region.End > region.Start, "End should be greater than start")
		assertEqual(t, region.End-region.Start, region.Size, "Size should match end-start")
		assertEqual(t, len(region.Coverage), region.Size, "Coverage array length should match size")
	}
}
