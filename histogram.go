package coverageplot

import (
	"fmt"
	log "github.com/sirupsen/logrus"
)

type Histogram struct {
	sum          float64 // running total of measurements
	Buckets      []int   // bucketed values for histogram
	BucketsCumm  []int
	bucketLabels []float64
	valueCount   int // number of values recorded for single value
	bucketCount  int
}

func (h *Histogram) init(bucketLabels []float64) {
	//h.bucketLabels = make([]float64, len(bucketLabels))
	//for i, v := range bucketLabels {
	//	h.bucketLabels[i] = float64(v)
	//}
	log.Debug(fmt.Sprintf("Histogram with labels: %v", bucketLabels))
	h.bucketLabels = bucketLabels
	h.bucketCount = len(bucketLabels)
	h.Buckets = make([]int, h.bucketCount)
	h.BucketsCumm = make([]int, h.bucketCount)
}

func (h *Histogram) Add(value float64) {
	h.sum += value
	bucketIndex := h.getBucket(value)
	h.valueCount++
	h.Buckets[bucketIndex]++
	for i := 0; i <= bucketIndex; i++ {
		h.BucketsCumm[i]++
	}
}

func (h *Histogram) AddAll(values []float64) {
	for _, v := range values {
		h.Add(v)
	}
}

func (h *Histogram) getBucket(i float64) int {
	for index, value := range h.bucketLabels {
		if i >= value {
			continue
		} else {
			if index == 0 {
				return 0
			} else {
				return index - 1
			}

		}
	}
	return h.bucketCount - 1
}
