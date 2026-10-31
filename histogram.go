package coverageplot

import (
	"log/slog"
)

// Histogram maintains discrete and cumulative bucket counts for observed values.
type Histogram struct {
	sum          float64
	Buckets      []int
	BucketsCumm  []int
	bucketLabels []float64
	valueCount   int
	bucketCount  int
}

// NewHistogram creates and initializes a Histogram with bucket thresholds.
func NewHistogram(bucketLabels []float64) *Histogram {
	h := &Histogram{}
	h.init(bucketLabels)
	return h
}

func (h *Histogram) init(bucketLabels []float64) {
	slog.Debug("Initializing histogram", "labels", bucketLabels)
	h.bucketLabels = bucketLabels
	h.bucketCount = len(bucketLabels)
	h.Buckets = make([]int, h.bucketCount)
	h.BucketsCumm = make([]int, h.bucketCount)
}

// AddCount records multiple occurrences of a float64 measurement into the histogram.
func (h *Histogram) AddCount(value float64, count int) {
	if count <= 0 {
		return
	}
	h.sum += value * float64(count)
	bucketIndex := h.getBucket(value)
	h.valueCount += count
	h.Buckets[bucketIndex] += count
	for i := 0; i <= bucketIndex; i++ {
		h.BucketsCumm[i] += count
	}
}

// Add records a single float64 measurement into the histogram.
func (h *Histogram) Add(value float64) {
	h.AddCount(value, 1)
}

// AddAll records multiple float64 measurements into the histogram.
func (h *Histogram) AddAll(values []float64) {
	for _, v := range values {
		h.Add(v)
	}
}

func (h *Histogram) getBucket(val float64) int {
	for idx, bound := range h.bucketLabels {
		if val < bound {
			if idx == 0 {
				return 0
			}
			return idx - 1
		}
	}
	return h.bucketCount - 1
}
