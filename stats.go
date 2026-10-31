package coverageplot

import (
	"errors"
	"math"
	"slices"
)

// Quartiles represents the first, second (median), and third quartiles.
type Quartiles struct {
	Q1 float64
	Q2 float64
	Q3 float64
}

// Min returns the minimum value in data.
func Min(data []float64) (float64, error) {
	if len(data) == 0 {
		return 0, errors.New("empty dataset")
	}
	m := data[0]
	for _, v := range data[1:] {
		if v < m {
			m = v
		}
	}
	return m, nil
}

// Max returns the maximum value in data.
func Max(data []float64) (float64, error) {
	if len(data) == 0 {
		return 0, errors.New("empty dataset")
	}
	m := data[0]
	for _, v := range data[1:] {
		if v > m {
			m = v
		}
	}
	return m, nil
}

// Mean returns the arithmetic mean of data.
func Mean(data []float64) (float64, error) {
	if len(data) == 0 {
		return 0, errors.New("empty dataset")
	}
	var sum float64
	for _, v := range data {
		sum += v
	}
	return sum / float64(len(data)), nil
}

// StandardDeviation returns the population standard deviation.
func StandardDeviation(data []float64) (float64, error) {
	if len(data) == 0 {
		return 0, errors.New("empty dataset")
	}
	mean, err := Mean(data)
	if err != nil {
		return 0, err
	}
	var sumSq float64
	for _, v := range data {
		diff := v - mean
		sumSq += diff * diff
	}
	return math.Sqrt(sumSq / float64(len(data))), nil
}

// Median returns the median of a pre-sorted slice.
func median(sorted []float64) float64 {
	n := len(sorted)
	if n == 0 {
		return 0
	}
	if n%2 == 1 {
		return sorted[n/2]
	}
	return (sorted[n/2-1] + sorted[n/2]) / 2.0
}

// Quartile computes Q1, Q2, and Q3 using Tukey's hinges.
func Quartile(data []float64) (Quartiles, error) {
	if len(data) < 2 {
		return Quartiles{}, errors.New("dataset must contain at least 2 elements")
	}
	sorted := slices.Clone(data)
	slices.Sort(sorted)

	q2 := median(sorted)
	cut := len(sorted) / 2
	lower := sorted[:cut]
	var upper []float64
	if len(sorted)%2 == 0 {
		upper = sorted[cut:]
	} else {
		upper = sorted[cut+1:]
	}

	return Quartiles{
		Q1: median(lower),
		Q2: q2,
		Q3: median(upper),
	}, nil
}

