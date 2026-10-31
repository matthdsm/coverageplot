package coverageplot

import (
	"slices"
	"testing"
)

func TestHistogram(t *testing.T) {
	labels := []float64{0, 10, 20, 30, 40}
	h := NewHistogram(labels)

	if len(h.Buckets) != len(labels) {
		t.Fatalf("expected %d buckets, got %d", len(labels), len(h.Buckets))
	}

	h.Add(-5)
	h.Add(5)
	h.Add(10)
	h.Add(25)
	h.Add(30)
	h.Add(50)

	expectedBuckets := []int{2, 1, 1, 1, 1}
	if !slices.Equal(h.Buckets, expectedBuckets) {
		t.Fatalf("expected buckets %v, got %v", expectedBuckets, h.Buckets)
	}

	// Cumulative histogram: sum of buckets 0..i
	// index 0: 2
	// index 1: 2 + 1 = 3
	// index 2: 3 + 1 = 4
	// index 3: 4 + 1 = 5
	// index 4: 5 + 1 = 6
	expectedCumm := []int{2, 3, 4, 5, 6}
	if !slices.Equal(h.BucketsCumm, expectedCumm) {
		t.Fatalf("expected cumulative buckets %v, got %v", expectedCumm, h.BucketsCumm)
	}
}

func TestHistogramAddCount(t *testing.T) {
	labels := []float64{0, 10, 20}
	h1 := NewHistogram(labels)
	for i := 0; i < 20; i++ {
		h1.Add(15.0)
	}

	h2 := NewHistogram(labels)
	h2.AddCount(15.0, 20)

	if !slices.Equal(h1.Buckets, h2.Buckets) {
		t.Fatalf("buckets mismatch: h1=%v, h2=%v", h1.Buckets, h2.Buckets)
	}
	if !slices.Equal(h1.BucketsCumm, h2.BucketsCumm) {
		t.Fatalf("cumulative mismatch: h1=%v, h2=%v", h1.BucketsCumm, h2.BucketsCumm)
	}
}

func TestHistogramAddAll(t *testing.T) {
	labels := []float64{0, 10, 20}
	h := NewHistogram(labels)
	h.AddAll([]float64{5, 15, 25})

	if h.valueCount != 3 {
		t.Fatalf("expected valueCount 3, got %d", h.valueCount)
	}
	expected := []int{1, 1, 1}
	if !slices.Equal(h.Buckets, expected) {
		t.Fatalf("expected %v, got %v", expected, h.Buckets)
	}
	expectedCumm := []int{1, 2, 3}
	if !slices.Equal(h.BucketsCumm, expectedCumm) {
		t.Fatalf("expected cumulative %v, got %v", expectedCumm, h.BucketsCumm)
	}
}
