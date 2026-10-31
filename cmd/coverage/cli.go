package main

import (
	"errors"
	"fmt"
	"io"
	"strings"
)

// ExitCode represents POSIX CLI exit statuses.
type ExitCode int

const (
	ExitSuccess ExitCode = 0
	ExitFailure ExitCode = 1
	ExitUsage   ExitCode = 2
)

// StatisticsConfig holds parsed options for coverage calculation.
type StatisticsConfig struct {
	CoverageFile  string
	IntervalsFile string
	PanelName     string
	TotalOnly     bool
	Debug         bool
}

// CLI holds the top-level CLI execution state.
type CLI struct {
	ShowHelp    bool
	ShowVersion bool
	Config      *StatisticsConfig
}

// ParseCLI parses arguments according to POSIX and GNU flag conventions.
func ParseCLI(args []string) (*CLI, ExitCode, error) {
	if len(args) == 0 {
		return &CLI{ShowHelp: true}, ExitUsage, errors.New("missing required arguments")
	}

	cli := &CLI{}

	// Skip optional "statistics" subcommand token if provided
	if args[0] == "statistics" {
		args = args[1:]
		if len(args) == 0 {
			return &CLI{ShowHelp: true}, ExitUsage, errors.New("missing required arguments")
		}
	}

	cfg := &StatisticsConfig{
		PanelName: "NA",
	}

	for i := 0; i < len(args); i++ {
		arg := args[i]
		if arg == "-h" || arg == "--help" || arg == "help" {
			cli.ShowHelp = true
			return cli, ExitSuccess, nil
		}
		if arg == "-v" || arg == "--version" || arg == "version" {
			cli.ShowVersion = true
			return cli, ExitSuccess, nil
		}
		if arg == "-t" || arg == "--total" {
			cfg.TotalOnly = true
			continue
		}
		if arg == "-d" || arg == "--debug" {
			cfg.Debug = true
			continue
		}

		if strings.HasPrefix(arg, "-") && !strings.HasPrefix(arg, "--") && len(arg) > 2 {
			bundled := arg[1:]
			allFlags := true
			for _, r := range bundled {
				if r == 't' {
					cfg.TotalOnly = true
				} else if r == 'd' {
					cfg.Debug = true
				} else if r == 'h' {
					cli.ShowHelp = true
					return cli, ExitSuccess, nil
				} else if r == 'v' {
					cli.ShowVersion = true
					return cli, ExitSuccess, nil
				} else {
					allFlags = false
					break
				}
			}
			if allFlags {
				continue
			}
		}

		if val, nextI, ok, err := extractFlagValue(args, i, "-c", "--coverage"); ok {
			if err != nil {
				return nil, ExitUsage, err
			}
			cfg.CoverageFile = val
			i = nextI
			continue
		}

		if val, nextI, ok, err := extractFlagValue(args, i, "-i", "--intervals"); ok {
			if err != nil {
				return nil, ExitUsage, err
			}
			cfg.IntervalsFile = val
			i = nextI
			continue
		}

		if val, nextI, ok, err := extractFlagValue(args, i, "-p", "--panel"); ok {
			if err != nil {
				return nil, ExitUsage, err
			}
			cfg.PanelName = val
			i = nextI
			continue
		}

		return nil, ExitUsage, fmt.Errorf("unrecognized option %q", arg)
	}

	if cfg.CoverageFile == "" {
		return nil, ExitUsage, errors.New("missing required option: -c, --coverage <file>")
	}
	if !cfg.TotalOnly && cfg.IntervalsFile == "" {
		return nil, ExitUsage, errors.New("missing required option: -i, --intervals <file> (required unless -t, --total)")
	}

	cli.Config = cfg
	return cli, ExitSuccess, nil
}

func extractFlagValue(args []string, i int, shortFlag, longFlag string) (string, int, bool, error) {
	arg := args[i]

	if arg == shortFlag || arg == longFlag {
		if i+1 >= len(args) || strings.HasPrefix(args[i+1], "-") {
			return "", i, true, fmt.Errorf("option requires an argument: %s", arg)
		}
		return args[i+1], i + 1, true, nil
	}

	if strings.HasPrefix(arg, longFlag+"=") {
		return strings.TrimPrefix(arg, longFlag+"="), i, true, nil
	}

	if strings.HasPrefix(arg, shortFlag) && len(arg) > len(shortFlag) {
		val := strings.TrimPrefix(arg, shortFlag)
		val = strings.TrimPrefix(val, "=")
		return val, i, true, nil
	}

	return "", i, false, nil
}

// PrintUsage outputs standardized UNIX manual synopsis and help text.
func PrintUsage(w io.Writer) {
	fmt.Fprintln(w, "usage: coverage [statistics] -c <coverage.bed.gz> [-i <intervals.bed>] [options]")
	fmt.Fprintln(w, "")
	fmt.Fprintln(w, "Calculate regional and summary coverage statistics from BED files.")
	fmt.Fprintln(w, "")
	fmt.Fprintln(w, "options:")
	fmt.Fprintln(w, "  -c, --coverage <file>   Input coverage BED or bgzip BED file (indexed with .tbi/.csi) [Required]")
	fmt.Fprintln(w, "  -i, --intervals <file>  Input panel BED file [Required unless -t]")
	fmt.Fprintln(w, "  -p, --panel <name>      Panel name annotation (default: NA)")
	fmt.Fprintln(w, "  -t, --total             Calculate only whole-file total statistics")
	fmt.Fprintln(w, "  -d, --debug             Enable verbose debug logging")
	fmt.Fprintln(w, "  -v, --version           Print version information and exit")
	fmt.Fprintln(w, "  -h, --help              Show this help message")
}
