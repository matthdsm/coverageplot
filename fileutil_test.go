package coverageplot

import (
	"compress/gzip"
	"io"
	"os"
	"testing"
)

func TestOpenReaderPlain(t *testing.T) {
	tmp, err := os.CreateTemp("", "plain_*.bed")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmp.Name())

	content := "chr1\t100\t200\tgeneA\nchr2\t300\t400\tgeneB\n"
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatal(err)
	}
	tmp.Close()

	r, err := openReader(tmp.Name())
	if err != nil {
		t.Fatalf("openReader plain failed: %v", err)
	}
	defer r.Close()

	data, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll failed: %v", err)
	}
	if string(data) != content {
		t.Fatalf("expected %q, got %q", content, string(data))
	}
}

func TestOpenReaderGzip(t *testing.T) {
	tmp, err := os.CreateTemp("", "compressed_*.bed.gz")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmp.Name())

	content := "chr1\t100\t200\tgeneA\nchr2\t300\t400\tgeneB\n"
	gw := gzip.NewWriter(tmp)
	if _, err := gw.Write([]byte(content)); err != nil {
		t.Fatal(err)
	}
	gw.Close()
	tmp.Close()

	r, err := openReader(tmp.Name())
	if err != nil {
		t.Fatalf("openReader gzip failed: %v", err)
	}
	defer r.Close()

	data, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll failed: %v", err)
	}
	if string(data) != content {
		t.Fatalf("expected %q, got %q", content, string(data))
	}
}

func TestOpenReaderStdin(t *testing.T) {
	r, err := openReader("-")
	if err != nil {
		t.Fatalf("openReader stdin failed: %v", err)
	}
	if r == nil {
		t.Fatal("expected non-nil reader for stdin")
	}
}

func TestOpenReaderNonExistent(t *testing.T) {
	_, err := openReader("/non/existent/path/for/sure.bed")
	if err == nil {
		t.Fatal("expected error opening non-existent file, got nil")
	}
}

func TestOpenReaderEmptyPath(t *testing.T) {
	r, err := openReader("")
	if err != nil {
		t.Fatalf("openReader empty string failed: %v", err)
	}
	if r == nil {
		t.Fatal("expected non-nil reader for empty path")
	}
}
