# coverageplot

`coverageplot` is a high-performance command-line utility for calculating regional and summary coverage statistics from genomic sequencing data formatted in BED or bgzip-compressed BED format.

Targeted at modern genomics pipelines, `coverageplot` operates with **zero external dependencies** using 100% Go standard library components, supports POSIX-compliant command-line flags, and integrates with `tabix` for fast indexed random-access queries.

---

## Features

- **Zero External Dependencies**: Built entirely with Go standard library (`log/slog`, `compress/gzip`, `slices`, `math`). No third-party modules or framework overhead.
- **Fast Indexed Queries**: Leverages `tabix` to query regions directly from bgzip-compressed BED files (`.bed.gz` with `.tbi` or `.csi` indexes).
- **POSIX-Compliant CLI**: Supports standard short flags (`-c`, `-i`, `-p`, `-t`, `-d`, `-v`, `-h`), GNU long options (`--coverage`, `--intervals`, etc.), bundled flags (e.g. `-td`), and standard exit codes (`0` success, `1` runtime error, `2` syntax error).
- **Numerically Stable One-Pass Statistics**:
  - Incremental Knuth / Chan-Golub-LeVeque $O(1)$ batch variance updates for running means and standard deviations.
  - Exact-frequency stream quantiles ($Q_1, Q_2, Q_3$).
  - Tukey's hinges quartiles and bucketed cumulative histograms.
- **Structured JSON Output**: Emits structured JSON to stdout suitable for downstream ingestion, reporting, or visualisation pipelines.

---

## Prerequisites

- **Go**: Version 1.21 or later (tested against Go 1.27.3).
- **HTSlib / tabix**: Required on `$PATH` for indexed regional queries (`tabix -R <intervals> <coverage.bed.gz>`).

---

## Installation

Clone the repository and build the binary:

```bash
git clone https://github.com/matthdsm/coverageplot.git
cd coverageplot

# Build local binary
go build -o bin/coverage ./cmd/coverage

# Or install to $GOPATH/bin
go install ./cmd/coverage
```

---

## Command-Line Usage

```text
usage: coverage [statistics] -c <coverage.bed.gz> [-i <intervals.bed>] [options]

Calculate regional and summary coverage statistics from BED files.

options:
  -c, --coverage <file>   Input coverage BED or bgzip BED file (indexed with .tbi/.csi) [Required]
  -i, --intervals <file>  Input panel BED file [Required unless -t]
  -p, --panel <name>      Panel name annotation (default: NA)
  -t, --total             Calculate only whole-file total statistics
  -d, --debug             Enable verbose debug logging
  -v, --version           Print version information and exit
  -h, --help              Show this help message
```

### Exit Codes
- `0`: Success.
- `1`: Processing or execution failure (e.g., file not found, tabix process error).
- `2`: Command-line syntax or usage error (e.g., missing required flag, invalid argument).

---

## Examples

### 1. Regional Statistics for a Gene Panel
Calculate per-exon and total coverage metrics for a targeted gene panel against a tabix-indexed coverage file:

```bash
coverage -c sample.coverage.bed.gz -i target_panel.bed -p HereditaryCancer > coverage_stats.json
```

### 2. Whole-File Total Statistics
Calculate summary statistics across an entire coverage BED file without regional decomposition:

```bash
coverage -c sample.coverage.bed.gz -t -p WholeGenome > total_stats.json
```

### 3. Verbose Debug Logging
Enable debug level logging output to stderr:

```bash
coverage -c sample.coverage.bed.gz -i target_panel.bed -d
```

---

## Input Formats

### Coverage File (`-c, --coverage`)
- Plain text BED or bgzip-compressed BED (`.bed.gz`).
- Format: `chrom`, `start`, `end`, `coverage_depth` (tab-delimited, 0-based half-open coordinates).
- For regional calculations, the file must be indexed with `tabix -p bed sample.coverage.bed.gz`.

### Intervals File (`-i, --intervals`)
- Standard 0-based half-open BED format containing target regions:
  ```text
  chr1    95      125     BRCA1;NM_007294.4;ENST00000357654;ENSE00001815124;1
  ```
- The 4th column optionally contains pipeline metadata formatted as `Gene;RefSeq;EnsemblTranscript;EnsemblExon;Rank`.

---

## JSON Output Structure

The output emitted to stdout matches the following schema:

```json
{
  "Name": "HereditaryCancer",
  "Total": {
    "Gene_panel": "HereditaryCancer",
    "Genome_build": "GRCh38",
    "Start": 0,
    "End": 90,
    "Size": 90,
    "Min": 10.0,
    "Max": 55.0,
    "Mean": 34.55,
    "StdDev": 14.22,
    "Quartiles": {
      "Q1": 25.0,
      "Q2": 35.0,
      "Q3": 50.0
    },
    "Zero_bases": 0,
    "Histogram": [0, 0, 15, 25, 30, 20, 0, 0, 0, 0, 0, 0],
    "CummulativeHistogram": [90, 90, 75, 50, 20, 0, 0, 0, 0, 0, 0, 0]
  },
  "Regions": [
    {
      "Gene_panel": "HereditaryCancer",
      "Genome_build": "GRCh38",
      "Chromosome": "1",
      "Start": 95,
      "End": 125,
      "Gene": "BRCA1",
      "Transcript_ensembl": "ENST00000357654",
      "Transcript_refseq": "NM_007294.4",
      "Exon_ensemble": "ENSE00001815124",
      "Rank": "1",
      "Size": 30,
      "Coverage": [10, 10, ...],
      "Min": 10.0,
      "Max": 50.0,
      "Mean": 33.33,
      "StdDev": 12.47,
      "Quartiles": { "Q1": 30.0, "Q2": 35.0, "Q3": 40.0 },
      "Zero_bases": 0,
      "Histogram": [...],
      "CummulativeHistogram": [...]
    }
  ]
}
```

---

## Testing

Run the full automated test suite:

```bash
# Run unit tests
go test -v ./...

# Run tests with race detection
go test -race ./...

# Run static analysis
go vet ./...
```

---

## License

This project is licensed under the MIT License or terms defined in repository root.
