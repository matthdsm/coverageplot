package coverageplot

import (
	"math"
	"testing"
)

func TestStatsCalculations(t *testing.T) {
	data := []float64{10, 20, 30, 40, 50}

	minVal, err := Min(data)
	if err != nil || minVal != 10 {
		t.Fatalf("Min: expected 10, got %v (err: %v)", minVal, err)
	}

	maxVal, err := Max(data)
	if err != nil || maxVal != 50 {
		t.Fatalf("Max: expected 50, got %v (err: %v)", maxVal, err)
	}

	meanVal, err := Mean(data)
	if err != nil || meanVal != 30 {
		t.Fatalf("Mean: expected 30, got %v (err: %v)", meanVal, err)
	}

	stdVal, err := StandardDeviation(data)
	if err != nil {
		t.Fatalf("StandardDeviation error: %v", err)
	}
	expectedStd := math.Sqrt(200.0) // population variance: sum((x-30)^2)/5 = (400+100+0+100+400)/5 = 1000/5 = 200
	if math.Abs(stdVal-expectedStd) > 1e-6 {
		t.Fatalf("StandardDeviation: expected %v, got %v", expectedStd, stdVal)
	}

	// Quartiles for 5 items: [10, 20, 30, 40, 50]
	// median = 30
	// lower = [10, 20], median of lower = 15
	// upper = [40, 50], median of upper = 45
	q, err := Quartile(data)
	if err != nil {
		t.Fatalf("Quartile error: %v", err)
	}
	if q.Q1 != 15 || q.Q2 != 30 || q.Q3 != 45 {
		t.Fatalf("Quartile: expected Q1=15, Q2=30, Q3=45, got %+v", q)
	}
}

func TestStatsEmpty(t *testing.T) {
	empty := []float64{}
	if _, err := Min(empty); err == nil {
		t.Fatal("expected error on Min(empty)")
	}
	if _, err := Max(empty); err == nil {
		t.Fatal("expected error on Max(empty)")
	}
	if _, err := Mean(empty); err == nil {
		t.Fatal("expected error on Mean(empty)")
	}
	if _, err := StandardDeviation(empty); err == nil {
		t.Fatal("expected error on StandardDeviation(empty)")
	}
	if _, err := Quartile(empty); err == nil {
		t.Fatal("expected error on Quartile(empty)")
	}
}

