package coverageplot

import (
	"flag"
	"fmt"
	"io/ioutil"
	"os"
	"os/exec"
	"strings"
	"testing"

	log "github.com/sirupsen/logrus"
	"github.com/stretchr/testify/assert"
)

// Global flag for debug mode
var debugTest bool

// TestMain handles custom initialization for tests
func TestMain(m *testing.M) {
	// Define and parse the debug flag
	flag.BoolVar(&debugTest, "debug", false, "Enable debug logging for tests")
	flag.Parse()

	// Set log level based on flag
	if debugTest {
		log.SetLevel(log.DebugLevel)
		log.Debug("Debug logging enabled for tests")
	} else {
		log.SetLevel(log.InfoLevel)
	}

	// Run the tests
	os.Exit(m.Run())
}

// Mock for exec.Command
func mockExecCommand(command string, args ...string) *exec.Cmd {
	cs := []string{"-test.run=TestHelperProcess", "--", command}
	cs = append(cs, args...)
	cmd := exec.Command(os.Args[0], cs...)

	// Add custom identifier for the test
	if len(args) > 1 && strings.Contains(args[1], "identical_end_test.bed.gz") {
		cmd.Env = []string{"GO_WANT_HELPER_PROCESS=1", "TEST_TYPE=identical_end"}
	} else if len(args) > 1 && strings.Contains(args[1], "overlap_test.bed.gz") {
		cmd.Env = []string{"GO_WANT_HELPER_PROCESS=1", "TEST_TYPE=overlap"}
	} else if len(args) > 1 && strings.Contains(args[1], "loop_termination_test.bed.gz") {
		cmd.Env = []string{"GO_WANT_HELPER_PROCESS=1", "TEST_TYPE=loop_termination"}
	} else {
		cmd.Env = []string{"GO_WANT_HELPER_PROCESS=1", "TEST_TYPE=standard"}
	}

	return cmd
}

// TestHelperProcess isn't a real test - it's used to mock exec.Command
func TestHelperProcess(t *testing.T) {
	if os.Getenv("GO_WANT_HELPER_PROCESS") != "1" {
		return
	}

	// Select mock output based on test type
	testType := os.Getenv("TEST_TYPE")

	if testType == "identical_end" {
		// Custom mock data for identical end test with different coverage values
		mockOutput := `chr1	90	100	15
chr1	100	110	35
chr1	110	130	30
chr1	130	150	55
chr1	110	130	30
chr1	130	150	55
`
		os.Stdout.WriteString(mockOutput)
	} else if testType == "overlap" {
		// Custom mock data for overlap test
		mockOutput := `chr1	90	100	10
chr1	100	110	30
chr1	110	120	40
chr1	120	150	50
chr1	120	150	50
`
		os.Stdout.WriteString(mockOutput)
	} else if testType == "real_data" {
		// Try to read from the TabixOutput_exome.bed file
		data, err := ioutil.ReadFile("data/TabixOutput_exome.bed")
		if err == nil {
			os.Stdout.Write(data)
		}
	} else if testType == "loop_termination" {
		// Custom mock data with critical segments at interval boundaries
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
	} else if testType == "tango2_data" {
		// Try to read from the TANGO2 TabixOutput file
		data, err := ioutil.ReadFile("data/tabixOut_TANGO2_test.exome.bed")
		if err == nil {
			os.Stdout.Write(data)
		} else {
			// Fallback mock data if file not found
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
	} else {
		// Standard mock output for other tests
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

/* Coordinate system conversion explanation:
 *
 * Tabix (1-based inclusive):  [ 101   102   103   104   105   ]
 * BedInterval (0-based):      [ 100   101   102   103   104   105)
 * Array indices:              [  0     1     2     3     4      ]
 *
 * Calculation: idx = i - interval.start
 *              where i is the 0-based position (start0 <= i < end0)
 */

func TestCalculateStatistics(t *testing.T) {
	// Save original exec.Command and restore after test
	origExecCommand := execCommand
	execCommand = mockExecCommand
	defer func() { execCommand = origExecCommand }()

	// Create a temporary bed file for testing
	tmpfile, err := ioutil.TempFile("", "test.bed")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmpfile.Name())

	// Setup test data - Fixed all Rank fields to be strings
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
			info:  []BedInfo{{Gene: "GENE2", TranscriptId: "ENST02", RefSeqId: "NM_002", ExonId: "ENSE02", Rank: "1"}}, // Changed 1 to "1"
		},
		{
			chrom: "chr2",
			start: 100,
			end:   130,
			info:  []BedInfo{{Gene: "GENE3", TranscriptId: "ENST03", RefSeqId: "NM_003", ExonId: "ENSE03", Rank: "1"}}, // Changed 1 to "1"
		},
	}

	// Run the function
	result := CalculateStatistics(intervals, tmpfile.Name(), "test.bed.gz", "TestPanel", false)
	log.Info(fmt.Sprintf("Coverage for region 1 is %v", result.Regions[1].Coverage))
	log.Info(fmt.Sprintf("Coverage for region 2 is %v", result.Regions[2].Coverage))
	// Verify results
	assert.Equal(t, []float64{10, 10, 10, 10, 10, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 50, 50, 50, 50, 50}, result.Regions[0].Coverage)
	assert.Equal(t, "TestPanel", result.Name)
	assert.Equal(t, 3, len(result.Regions)) // Check first interval statistics
	assert.Equal(t, "GENE1", result.Regions[0].Gene)
	assert.Equal(t, 95, result.Regions[0].Start)
	assert.Equal(t, 125, result.Regions[0].End)
	assert.Equal(t, 30, result.Regions[0].Size)
	assert.Equal(t, 10.0, result.Regions[0].Min)
	assert.Equal(t, 50.0, result.Regions[0].Max)
	assert.InEpsilon(t, 33.33, result.Regions[0].Mean, 0.01)

	// Check second interval statistics
	assert.Equal(t, []float64{20, 20, 20, 20, 20, 20, 20, 20, 20, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30}, result.Regions[1].Coverage)
	assert.Equal(t, "GENE2", result.Regions[1].Gene)
	assert.Equal(t, 201, result.Regions[1].Start)
	assert.Equal(t, 231, result.Regions[1].End)
	assert.Equal(t, 30, result.Regions[1].Size)
	assert.Equal(t, 20.0, result.Regions[1].Min)
	assert.Equal(t, 30.0, result.Regions[1].Max)
	assert.InEpsilon(t, 25.33, result.Regions[1].Mean, 0.01)

	// Check third interval statistics
	assert.Equal(t, []float64{45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 0, 0, 0, 50, 50, 50, 50, 50, 50, 50, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55}, result.Regions[2].Coverage)
	assert.Equal(t, "GENE3", result.Regions[2].Gene)
	assert.Equal(t, 100, result.Regions[2].Start)
	assert.Equal(t, 130, result.Regions[2].End)
	assert.Equal(t, 30, result.Regions[2].Size)
	assert.Equal(t, 0.0, result.Regions[2].Min)
	assert.Equal(t, 55.0, result.Regions[2].Max)
	assert.InEpsilon(t, 45.0, result.Regions[2].Mean, 0.01)

	// Check total statistics
	//log.Info(fmt.Sprintf("Coverage for total is %v", result.Total.Coverage))
	assert.Equal(t, 90, result.Total.Size) // 3 intervals of 30 bases each
	assert.Equal(t, 0.0, result.Total.Min)
	assert.Equal(t, 55.0, result.Total.Max)

	// Test with total flag
	totalResult := CalculateStatistics(intervals, tmpfile.Name(), "test.bed.gz", "TestPanel", true)
	assert.Equal(t, 0, len(totalResult.Regions))
	assert.Equal(t, 90, totalResult.Total.Size)
}

func TestCalculateStatisticsWithOverlappingIntervals(t *testing.T) {
	// Save original exec.Command and restore after test
	origExecCommand := execCommand
	execCommand = mockExecCommand
	defer func() { execCommand = origExecCommand }()

	// Create a temporary bed file for testing
	tmpfile, err := ioutil.TempFile("", "test.bed")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmpfile.Name())

	// Setup test data with overlapping intervals - Fixed all Rank fields to be strings
	intervals := []BedInterval{
		{
			chrom: "chr1",
			start: 95,
			end:   125,
			info:  []BedInfo{{Gene: "GENE1", TranscriptId: "ENST01", RefSeqId: "NM_001", ExonId: "ENSE01", Rank: "1"}}, // Changed 1 to "1"
		},
		{
			chrom: "chr1",
			start: 120,
			end:   150,
			info:  []BedInfo{{Gene: "GENE2", TranscriptId: "ENST02", RefSeqId: "NM_002", ExonId: "ENSE02", Rank: "1"}}, // Changed 1 to "1"
		},
	}

	// Run the function
	result := CalculateStatistics(intervals, tmpfile.Name(), "overlap_test.bed.gz", "TestPanel", false)

	// Verify results account for overlapping regions correctly
	assert.Equal(t, 2, len(result.Regions))
	assert.Equal(t, 55, result.Total.Size) // Should be 55 bases total (not 60) due to overlap
	assert.Equal(t, []float64{10, 10, 10, 10, 10, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 50, 50, 50, 50, 50}, result.Regions[0].Coverage)
	assert.Equal(t, 50.0, result.Regions[1].Min)
	assert.Equal(t, 50.0, result.Regions[1].Max)
}

func TestCalculateStatisticsWithOverlappingIdenticalEndIntervals(t *testing.T) {
	// Save original exec.Command and restore after test
	origExecCommand := execCommand
	execCommand = mockExecCommand
	defer func() { execCommand = origExecCommand }()

	// Create a temporary bed file for testing
	tmpfile, err := ioutil.TempFile("", "test.bed")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmpfile.Name())

	// Setup test data with overlapping intervals
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

	// Use a different filename to trigger the custom mock data
	result := CalculateStatistics(intervals, tmpfile.Name(), "identical_end_test.bed.gz", "TestPanel", false)

	// Check if second region is correct ( will be zeros if loop breaks at the wrong time
	assert.Equal(t, []float64{30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55}, result.Regions[1].Coverage)

	// Verify results account for overlapping regions correctly
	assert.Equal(t, 2, len(result.Regions))
	assert.Equal(t, 55, result.Total.Size)
	assert.Equal(t, 55.0, result.Total.Max)
	assert.Equal(t, 15.0, result.Total.Min)
}

func TestCalculateStatisticsWithEdgeCases(t *testing.T) {
	// Similar setup to previous tests

	// Test with empty intervals slice
	emptyResult := CalculateStatistics([]BedInterval{}, "intervals", "test.bed.gz", "EmptyPanel", false)
	assert.Equal(t, 0, len(emptyResult.Regions))
	assert.Equal(t, 0, emptyResult.Total.Size)

	// Additional edge cases could be tested here
}

func TestCalculateStatisticsWithRealData(t *testing.T) {
	// Save original exec.Command and restore after test
	origExecCommand := execCommand

	// Temporarily set log level to InfoLevel for this test
	origLogLevel := log.GetLevel()
	log.SetLevel(log.InfoLevel)
	defer log.SetLevel(origLogLevel)

	// Create a mock for this specific test that returns content from TabixOutput_exome.bed
	execCommand = func(command string, args ...string) *exec.Cmd {
		cs := []string{"-test.run=TestHelperProcess", "--", command}
		cs = append(cs, args...)
		cmd := exec.Command(os.Args[0], cs...)
		cmd.Env = []string{"GO_WANT_HELPER_PROCESS=1", "TEST_TYPE=real_data"}
		return cmd
	}
	defer func() { execCommand = origExecCommand }()

	// Set up intervals file or use ParseRegions
	intervalsFile := "data/intervals_per_exon.bed"
	var intervals []BedInterval

	// Check if the real intervals file exists
	if _, err := os.Stat(intervalsFile); os.IsNotExist(err) {
		t.Skip("Skipping test: intervals file not found")
		return
	}

	// Parse real intervals
	intervals = ParseRegions(intervalsFile)

	// Verify we got some intervals
	if len(intervals) == 0 {
		t.Fatal("No intervals parsed from the real intervals file")
	}

	// Log the first few intervals for debugging
	for i, intv := range intervals[:5] {
		t.Logf("Interval %d: %s:%d-%d, Gene: %s",
			i, intv.chrom, intv.start, intv.end,
			intv.info[0].Gene)
	}

	// Calculate statistics using our mocked tabix output
	result := CalculateStatistics(intervals, intervalsFile, "real_data_test.bed.gz", "RealDataTest", false)

	// Basic verification
	assert.Equal(t, "RealDataTest", result.Name)
	assert.True(t, len(result.Regions) > 0, "Should have at least one region")
	assert.True(t, result.Total.Size > 0, "Total size should be greater than 0")
	assert.Equal(t, 5.0, result.Total.Min)
	assert.Equal(t, 259.0, result.Total.Max)
	assert.InEpsilon(t, 97.72, result.Total.Mean, 0.01)

	// Log some results for manual verification
	t.Logf("Total regions: %d", len(result.Regions))
	t.Logf("Total covered bases: %d", result.Total.Size)
	t.Logf("Min coverage: %.2f", result.Total.Min)
	t.Logf("Max coverage: %.2f", result.Total.Max)
	t.Logf("Mean coverage: %.2f", result.Total.Mean)

	// Check that regions were properly processed
	for i, region := range result.Regions[:5] { // Check first 5 regions
		t.Logf("Region %d: %s:%d-%d, Min: %.2f, Max: %.2f, Mean: %.2f",
			i, region.Chromosome, region.Start, region.End,
			region.Min, region.Max, region.Mean)

		// These should be true for valid regions
		assert.True(t, region.End > region.Start, "End should be greater than start")
		assert.Equal(t, region.End-region.Start, region.Size, "Size should match end-start")
		assert.Equal(t, len(region.Coverage), region.Size, "Coverage array length should match size")
	}
}

func TestCalculateStatisticsWithCriticalLoopTermination(t *testing.T) {
	// Save original exec.Command and restore after test
	origExecCommand := execCommand

	// Create a specialized mock that tests the loop termination condition
	execCommand = func(command string, args ...string) *exec.Cmd {
		cs := []string{"-test.run=TestHelperProcess", "--", command}
		cs = append(cs, args...)
		cmd := exec.Command(os.Args[0], cs...)
		cmd.Env = []string{"GO_WANT_HELPER_PROCESS=1", "TEST_TYPE=loop_termination"}
		return cmd
	}
	defer func() { execCommand = origExecCommand }()

	// Create a temporary bed file for testing
	tmpfile, err := ioutil.TempFile("", "test.bed")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmpfile.Name())

	// Setup intervals to test critical boundary conditions
	intervals := []BedInterval{
		{
			// Test 1: Where a segment ends exactly at interval end (segment.end == interval.end)
			chrom: "chr1",
			start: 95,
			end:   125,
			info:  []BedInfo{{Gene: "BOUNDARY1", TranscriptId: "ENST01", RefSeqId: "NM_001", ExonId: "ENSE01", Rank: "1"}},
		},
		{
			// Test 2: Where we have zero-length segment and a segment that starts at boundary
			chrom: "chr2",
			start: 95,
			end:   125,
			info:  []BedInfo{{Gene: "BOUNDARY2", TranscriptId: "ENST02", RefSeqId: "NM_002", ExonId: "ENSE02", Rank: "1"}},
		},
	}

	// Run the function
	result := CalculateStatistics(intervals, tmpfile.Name(), "loop_termination.bed.gz", "TerminationTest", false)

	// TEST CASE 1: Check that coverage at the end boundary was correctly captured
	// This would fail if using "if end > interval.end" because the segment exactly ending at the boundary would be missed
	assert.Equal(t, 10.0, result.Regions[0].Coverage[0])  // First base should have coverage=10
	assert.Equal(t, 40.0, result.Regions[0].Coverage[29]) // Last base should have coverage=40
	assert.Equal(t, 40.0, result.Regions[0].Max)          // Max should be 40 (last segment)

	// TEST CASE 2: Check that zero-length segment at boundary was handled correctly
	// This would fail if using "if end > interval.end" because the 125-125 segment wouldn't trigger the break
	assert.Equal(t, 60.0, result.Regions[1].Coverage[29]) // Last base should have coverage=60, not 70 or 80

	// Test total aggregation is correct
	assert.Equal(t, 2, len(result.Regions))
	assert.Equal(t, 60, result.Total.Size)  // 30 + 30
	assert.Equal(t, 60.0, result.Total.Max) // Overall max
	assert.Equal(t, 10.0, result.Total.Min) // Overall min
}

func TestCalculateStatisticsWithTango2Data(t *testing.T) {
	// Save original exec.Command and restore after test
	origExecCommand := execCommand

	// Create a mock for this specific test that returns content from TANGO2 bed file
	execCommand = func(command string, args ...string) *exec.Cmd {
		cs := []string{"-test.run=TestHelperProcess", "--", command}
		cs = append(cs, args...)
		cmd := exec.Command(os.Args[0], cs...)
		cmd.Env = []string{"GO_WANT_HELPER_PROCESS=1", "TEST_TYPE=tango2_data"}
		return cmd
	}
	defer func() { execCommand = origExecCommand }()

	// Set up intervals file for TANGO2
	intervalsFile := "data/TANGO2_intervals_per_exon.bed"
	var intervals []BedInterval

	// Check if the TANGO2 intervals file exists
	if _, err := os.Stat(intervalsFile); os.IsNotExist(err) {
		t.Skip("Skipping test: TANGO2 intervals file not found")
		return
	}

	// Parse TANGO2 intervals
	intervals = ParseRegions(intervalsFile)

	// Verify we got some intervals
	if len(intervals) == 0 {
		t.Fatal("No intervals parsed from the TANGO2 intervals file")
	}

	// Log the first few intervals for debugging
	for i, intv := range intervals[3:6] {
		t.Logf("TANGO2 Interval %d: %s:%d-%d, Gene: %s",
			i, intv.chrom, intv.start, intv.end,
			intv.info[0].Gene)
	}

	// Calculate statistics using our mocked tabix output
	result := CalculateStatistics(intervals, intervalsFile, "data/tabixOut_TANGO2_test.exome.bed", "TANGO2Test", false)

	// Basic verification
	assert.Equal(t, "TANGO2Test", result.Name)
	assert.True(t, len(result.Regions) > 0, "Should have at least one region")
	assert.True(t, result.Total.Size > 0, "Total size should be greater than 0")

	// Log some results for manual verification
	t.Logf("TANGO2 Total regions: %d", len(result.Regions))
	t.Logf("TANGO2 Total covered bases: %d", result.Total.Size)
	t.Logf("TANGO2 Min coverage: %.2f", result.Total.Min)
	t.Logf("TANGO2 Max coverage: %.2f", result.Total.Max)
	t.Logf("TANGO2 Mean coverage: %.2f", result.Total.Mean)

	// Check region i=3
	assert.InEpsilon(t, 73.8, result.Regions[3].Mean, 0.01)
	assert.InEpsilon(t, 77.23, result.Regions[4].Mean, 0.01)
	assert.InEpsilon(t, 108.41, result.Regions[5].Mean, 0.01)

	// Check that regions were properly processed
	for i, region := range result.Regions[3:6] {
		t.Logf("TANGO2 Region %d: %s:%d-%d, Min: %.2f, Max: %.2f, Mean: %.2f",
			i, region.Chromosome, region.Start, region.End,
			region.Min, region.Max, region.Mean)

		// These should be true for valid regions
		assert.True(t, region.End > region.Start, "End should be greater than start")
		assert.Equal(t, region.End-region.Start, region.Size, "Size should match end-start")
		assert.Equal(t, len(region.Coverage), region.Size, "Coverage array length should match size")
	}
}
