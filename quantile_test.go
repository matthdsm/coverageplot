package coverageplot

import (
	"testing"
)

func TestStreamQuantiles(t *testing.T) {
	q := NewTargeted(0.25, 0.50, 0.75)

	for i := 1; i <= 100; i++ {
		q.Insert(float64(i))
	}

	q25 := q.Query(0.25)
	q50 := q.Query(0.50)
	q75 := q.Query(0.75)

	if q25 != 25.0 {
		t.Fatalf("expected Q25=25, got %v", q25)
	}
	if q50 != 50.0 {
		t.Fatalf("expected Q50=50, got %v", q50)
	}
	if q75 != 75.0 {
		t.Fatalf("expected Q75=75, got %v", q75)
	}
}

func TestStreamAddCount(t *testing.T) {
	q := NewTargeted(0.50)
	q.AddCount(10.0, 50)
	q.AddCount(20.0, 50)

	median := q.Query(0.50)
	if median != 10.0 {
		t.Fatalf("expected median 10.0, got %v", median)
	}
}

func TestStreamEmpty(t *testing.T) {
	q := NewTargeted(0.50)
	if v := q.Query(0.50); v != 0.0 {
		t.Fatalf("expected 0.0 for empty stream, got %v", v)
	}
}

