package main

import (
	"bytes"
	"strings"
	"testing"
)

func TestParseCLISuccess(t *testing.T) {
	// Direct execution with short flags
	cli, code, err := ParseCLI([]string{"-c", "cov.bed.gz", "-i", "panel.bed", "-p", "BRCA", "-t", "-d"})
	if err != nil || code != ExitSuccess {
		t.Fatalf("unexpected error: %v, code: %d", err, code)
	}
	if cli.Config.CoverageFile != "cov.bed.gz" || cli.Config.IntervalsFile != "panel.bed" || cli.Config.PanelName != "BRCA" || !cli.Config.TotalOnly || !cli.Config.Debug {
		t.Fatalf("unexpected config: %+v", cli.Config)
	}

	// Subcommand keyword with long flags and = syntax
	cli, code, err = ParseCLI([]string{"statistics", "--coverage=cov.bed.gz", "--total", "--debug"})
	if err != nil || code != ExitSuccess {
		t.Fatalf("unexpected error: %v, code: %d", err, code)
	}
	if cli.Config.CoverageFile != "cov.bed.gz" || !cli.Config.TotalOnly || !cli.Config.Debug {
		t.Fatalf("unexpected config: %+v", cli.Config)
	}
}

func TestParseCLIBundledFlags(t *testing.T) {
	cli, code, err := ParseCLI([]string{"-c", "cov.bed.gz", "-td"})
	if err != nil || code != ExitSuccess {
		t.Fatalf("unexpected error: %v, code: %d", err, code)
	}
	if !cli.Config.TotalOnly || !cli.Config.Debug {
		t.Fatalf("expected TotalOnly and Debug true for bundled -td, got %+v", cli.Config)
	}
}

func TestParseCLIErrors(t *testing.T) {
	// Empty args
	_, code, err := ParseCLI([]string{})
	if err == nil || code != ExitUsage {
		t.Fatalf("expected ExitUsage on empty args, got code: %d, err: %v", code, err)
	}

	// Missing coverage flag
	_, code, err = ParseCLI([]string{"-i", "panel.bed"})
	if err == nil || code != ExitUsage {
		t.Fatalf("expected ExitUsage on missing coverage flag, got code: %d, err: %v", code, err)
	}

	// Missing intervals when not total
	_, code, err = ParseCLI([]string{"-c", "cov.bed.gz"})
	if err == nil || code != ExitUsage {
		t.Fatalf("expected ExitUsage on missing intervals flag, got code: %d, err: %v", code, err)
	}

	// Flag requiring value missing value
	_, code, err = ParseCLI([]string{"-c"})
	if err == nil || code != ExitUsage {
		t.Fatalf("expected ExitUsage on missing flag argument, got code: %d, err: %v", code, err)
	}
}

func TestParseCLIHelpAndVersion(t *testing.T) {
	cli, code, err := ParseCLI([]string{"--help"})
	if err != nil || code != ExitSuccess || !cli.ShowHelp {
		t.Fatalf("expected ShowHelp true on --help, got: %+v, code: %d, err: %v", cli, code, err)
	}

	cli, code, err = ParseCLI([]string{"--version"})
	if err != nil || code != ExitSuccess || !cli.ShowVersion {
		t.Fatalf("expected ShowVersion true on --version, got: %+v, code: %d, err: %v", cli, code, err)
	}
}

func TestPrintUsage(t *testing.T) {
	var buf bytes.Buffer
	PrintUsage(&buf)
	out := buf.String()
	if !strings.Contains(out, "-c, --coverage") {
		t.Fatalf("usage missing expected flags: %s", out)
	}
}
