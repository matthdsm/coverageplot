package coverageplot

import (
	"slices"
)

// Stream computes quantiles over a data stream using frequency tracking.
type Stream struct {
	counts map[float64]uint64
	total  uint64
}

// NewTargeted returns a new quantile Stream for the requested target quantiles.
func NewTargeted(targets ...float64) *Stream {
	return &Stream{
		counts: make(map[float64]uint64),
	}
}

// Insert adds a value observation to the stream.
func (s *Stream) Insert(v float64) {
	s.counts[v]++
	s.total++
}

// AddCount records multiple occurrences of a value observation.
func (s *Stream) AddCount(v float64, n uint64) {
	if n == 0 {
		return
	}
	s.counts[v] += n
	s.total += n
}

// Query returns the value corresponding to quantile q (0.0 <= q <= 1.0).
func (s *Stream) Query(q float64) float64 {
	if s.total == 0 || len(s.counts) == 0 {
		return 0.0
	}
	targetRank := q * float64(s.total)

	keys := make([]float64, 0, len(s.counts))
	for k := range s.counts {
		keys = append(keys, k)
	}
	slices.Sort(keys)

	var cum uint64
	for _, val := range keys {
		cum += s.counts[val]
		if float64(cum) >= targetRank {
			return val
		}
	}
	return keys[len(keys)-1]
}

