package main

import (
	"fmt"
	"log/slog"
	"os"
	"strings"

	coverage "gitlab.cmgg.be/cmgg/coverageplot"
)

const version = "1.0.0"

func main() {
	cli, code, err := ParseCLI(os.Args[1:])
	if err != nil {
		fmt.Fprintf(os.Stderr, "coverage: %s\n\n", err)
		PrintUsage(os.Stderr)
		os.Exit(int(code))
	}

	if cli.ShowVersion {
		fmt.Printf("coverage version %s\n", version)
		os.Exit(int(ExitSuccess))
	}

	if cli.ShowHelp {
		PrintUsage(os.Stdout)
		os.Exit(int(ExitSuccess))
	}

	debugEnabled := false
	if env := os.Getenv("ENVIRONMENT"); strings.EqualFold(env, "development") {
		debugEnabled = true
	}
	if env := os.Getenv("LOG_LEVEL"); strings.EqualFold(env, "debug") {
		debugEnabled = true
	}
	if cli.Config != nil && cli.Config.Debug {
		debugEnabled = true
	}

	logLevel := slog.LevelInfo
	if debugEnabled {
		logLevel = slog.LevelDebug
	}
	handler := slog.NewTextHandler(os.Stderr, &slog.HandlerOptions{Level: logLevel})
	slog.SetDefault(slog.New(handler))

	err = coverage.Statistics(
		cli.Config.IntervalsFile,
		cli.Config.CoverageFile,
		cli.Config.PanelName,
		cli.Config.TotalOnly,
	)
	if err != nil {
		slog.Error("Execution failed", "err", err)
		os.Exit(int(ExitFailure))
	}
}
