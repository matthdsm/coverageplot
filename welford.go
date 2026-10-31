package coverageplot

import (
	"math"
)

// Stats gathers running statistics on values fed into it using Welford's algorithm.
type Stats struct {
	count             uint64
	min, max, mean, s float64
}

// New returns a new Stats instance.
func New() *Stats {
	return &Stats{}
}

// Min returns the minimum value seen.
func (s *Stats) Min() float64 {
	return s.min
}

// Max returns the maximum value seen.
func (s *Stats) Max() float64 {
	return s.max
}

// Count returns the total number of values seen.
func (s *Stats) Count() uint64 {
	return s.count
}

// Mean returns the mean of all values seen.
func (s *Stats) Mean() float64 {
	return s.mean
}

// Variance returns the population variance of all values seen.
func (s *Stats) Variance() float64 {
	if s.count > 0 {
		return s.s / float64(s.count)
	}
	return 0.0
}

// Stddev returns the population standard deviation of all values seen.
func (s *Stats) Stddev() float64 {
	return math.Sqrt(s.Variance())
}

// Reset sets all counters to zero.
func (s *Stats) Reset() {
	s.count = 0
	s.min = 0.0
	s.max = 0.0
	s.mean = 0.0
	s.s = 0.0
}

// AddCount feeds n occurrences of a value into Stats in O(1) time using Chan-Golub-LeVeque merging.
func (s *Stats) AddCount(val float64, n uint64) {
	if n == 0 {
		return
	}
	if s.count == 0 {
		s.min = val
		s.max = val
		s.mean = val
		s.s = 0.0
		s.count = n
		return
	}
	if val < s.min {
		s.min = val
	}
	if val > s.max {
		s.max = val
	}
	delta := val - s.mean
	newCount := s.count + n
	s.mean += delta * (float64(n) / float64(newCount))
	s.s += delta * delta * (float64(s.count) * float64(n) / float64(newCount))
	s.count = newCount
}

// Add feeds a single value into the Stats.
func (s *Stats) Add(val float64) {
	s.AddCount(val, 1)
}
