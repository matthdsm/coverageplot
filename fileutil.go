package coverageplot

import (
	"compress/gzip"
	"io"
	"os"
)

type gzipReadCloser struct {
	gz *gzip.Reader
	f  *os.File
}

func (g *gzipReadCloser) Read(p []byte) (int, error) {
	return g.gz.Read(p)
}

func (g *gzipReadCloser) Close() error {
	errGz := g.gz.Close()
	errF := g.f.Close()
	if errGz != nil {
		return errGz
	}
	return errF
}

// openReader returns a ReadCloser for the given file path.
// If path is "-" or empty, it returns os.Stdin.
// If the file is gzip/bgzf compressed, it transparently decompresses it.
func openReader(path string) (io.ReadCloser, error) {
	if path == "" || path == "-" {
		return io.NopCloser(os.Stdin), nil
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}

	var header [2]byte
	n, readErr := io.ReadFull(f, header[:])
	if readErr == nil && n == 2 && header[0] == 0x1f && header[1] == 0x8b {
		if _, seekErr := f.Seek(0, io.SeekStart); seekErr == nil {
			gz, gzErr := gzip.NewReader(f)
			if gzErr != nil {
				f.Close()
				return nil, gzErr
			}
			return &gzipReadCloser{gz: gz, f: f}, nil
		}
	}

	if _, seekErr := f.Seek(0, io.SeekStart); seekErr != nil {
		f.Close()
		return nil, seekErr
	}

	return f, nil
}
