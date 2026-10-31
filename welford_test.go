package coverageplot

import (
	"math"
	"testing"
)

func TestWelfordStats(t *testing.T) {
	s := New()
	values := []float64{2, 4, 4, 4, 5, 5, 7, 9}

	for _, v := range values {
		s.Add(v)
	}

	if s.Count() != 8 {
		t.Fatalf("expected count 8, got %d", s.Count())
	}
	if s.Min() != 2.0 {
		t.Fatalf("expected min 2.0, got %f", s.Min())
	}
	if s.Max() != 9.0 {
		t.Fatalf("expected max 9.0, got %f", s.Max())
	}
	if s.Mean() != 5.0 {
		t.Fatalf("expected mean 5.0, got %f", s.Mean())
	}
	expectedVar := 32.0 / 8.0
	if math.Abs(s.Variance()-expectedVar) > 1e-6 {
		t.Fatalf("expected variance %f, got %f", expectedVar, s.Variance())
	}
	if math.Abs(s.Stddev()-math.Sqrt(expectedVar)) > 1e-6 {
		t.Fatalf("expected stddev %f, got %f", math.Sqrt(expectedVar), s.Stddev())
	}

	s.Reset()
	if s.Count() != 0 || s.Mean() != 0.0 || s.Min() != 0.0 || s.Max() != 0.0 || s.Variance() != 0.0 {
		t.Fatalf("expected zeroed stats after Reset, got %+v", s)
	}
}

func TestWelfordBatchAddCount(t *testing.T) {
	s1 := New()
	for i := 0; i < 50; i++ {
		s1.Add(10.0)
	}
	for i := 0; i < 50; i++ {
		s1.Add(20.0)
	}

	s2 := New()
	s2.AddCount(10.0, 50)
	s2.AddCount(20.0, 50)

	if s1.Count() != s2.Count() {
		t.Fatalf("count mismatch: s1=%d, s2=%d", s1.Count(), s2.Count())
	}
	if math.Abs(s1.Mean()-s2.Mean()) > 1e-9 {
		t.Fatalf("mean mismatch: s1=%f, s2=%f", s1.Mean(), s2.Mean())
	}
	if math.Abs(s1.Variance()-s2.Variance()) > 1e-9 {
		t.Fatalf("variance mismatch: s1=%f, s2=%f", s1.Variance(), s2.Variance())
	}
}
