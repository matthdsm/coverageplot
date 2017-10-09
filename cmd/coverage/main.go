/*
The coverageplot command performs a coverage analysis on sequencing data.
You can use the d3.js based interface to visualise the results of the analysis
in an exon centric way.
*/
package main

import (
	"flag"
	"fmt"
	log "github.com/sirupsen/logrus"
	coverage "gitlab.cmgg.be/cmgg/coverageplot"
	"os"
	"strconv"
)

func main() {
	analyseCommand := flag.NewFlagSet("analyse", flag.ExitOnError)
	analyseIntervalsFlag := analyseCommand.String("intervals", "", "Input the sorted and merged interval file.")
	bamFlag := analyseCommand.String("bam", "", "Input bam file, a sorted and indexed file is required.")
	formatFlagAnalyse := analyseCommand.String("format", "bed", "format of the intervals file (bed|gff)")

	statisticsCommand := flag.NewFlagSet("statistics", flag.ExitOnError)
	statisticsIntervalsFlag := statisticsCommand.String("intervals", "", "Input the sorted panel bed file file.")
	statisticsCoverageFlag := statisticsCommand.String("coverage", "", "Input the sorted coverage bed file .")
	statisticsPanelFlag := statisticsCommand.String("panel", "", "Input a name for the panel (default: NA).")

	serveCommand := flag.NewFlagSet("serve", flag.ExitOnError)

	intervalsCommand := flag.NewFlagSet("intervals", flag.ExitOnError)
	gffFlag := intervalsCommand.String("gff", "", "Input gff file.")

	mergeCommand := flag.NewFlagSet("merge", flag.ExitOnError)
	intervalsFlag := mergeCommand.String("intervals", "", "Input the sorted interval file.")
	expandFlag := mergeCommand.String("expand", "", "Expand the interval with # bp on each side")
	formatFlag := mergeCommand.String("format", "bed", "format of the input and output files (bed|gff)")

	if len(os.Args) == 1 {
		fmt.Println(os.Stderr, "usage: coverage <command> [<args>]")
		fmt.Println(os.Stderr, "The available commands are: ")
		fmt.Println(os.Stderr, " analyse    Performs a coverage analysis")
		fmt.Println(os.Stderr, " serve      Run webserver loading the results")
		fmt.Println(os.Stderr, " intervals  Parse gff file and generate db of all exon intervals")
		fmt.Println(os.Stderr, " merge      Input the sorted interval file (sort -k1,1n -k2,2n)")
		return
	}

	switch os.Args[1] {
	case "analyse":
		analyseCommand.Parse(os.Args[2:])
	case "serve":
		serveCommand.Parse(os.Args[2:])
	case "intervals":
		intervalsCommand.Parse(os.Args[2:])
	case "merge":
		mergeCommand.Parse(os.Args[2:])
	case "statistics":
		statisticsCommand.Parse(os.Args[2:])
	default:
		log.Debug(fmt.Sprintf("%q is not valid command.\n", os.Args[1]))
		os.Exit(2)
	}

	Environment := os.Getenv("Environment")
	if Environment == "development" {
		log.SetLevel(log.DebugLevel)
	} else {
		// The TextFormatter is default, you don't actually have to do this.
		log.SetLevel(log.InfoLevel)
	}

	if statisticsCommand.Parsed() {
		if *statisticsIntervalsFlag == "" {
			fmt.Println(os.Stderr, "Please supply a sorted panel bed file using -intervals option.")
			os.Exit(2)
		}
		if *statisticsCoverageFlag == "" {
			fmt.Println(os.Stderr, "Please supply a sorted coverage bed file using -coverage option.")
			os.Exit(2)
		}
		if *statisticsPanelFlag == "" {
			*statisticsPanelFlag = "NA"
		}
		log.Info(fmt.Sprintf("Analysing file %q using intervals %q", *statisticsCoverageFlag, *statisticsIntervalsFlag))
		coverage.Statistics(*statisticsIntervalsFlag, *statisticsCoverageFlag, *statisticsPanelFlag)
	}

	if analyseCommand.Parsed() {
		if *bamFlag == "" {
			fmt.Println(os.Stderr, "Please supply a bam file as input using -bam option.")
			return
		}
		if *analyseIntervalsFlag == "" {
			fmt.Println(os.Stderr, "Please supply a sorted and merged interval file using -intervals option.")
			return
		}
		if *formatFlagAnalyse != "bed" && *formatFlagAnalyse != "gff" {
			fmt.Println(os.Stderr, "Please supply a valid format (bed|gff) using -format option.")
			os.Exit(2)
		}
		log.Info(fmt.Sprintf("Analysing bam file %q using intervals %q", *bamFlag, *analyseIntervalsFlag))
		coverage.Analyse(*bamFlag, *analyseIntervalsFlag, *formatFlagAnalyse)
	}

	if serveCommand.Parsed() {
		serve()
	}

	if intervalsCommand.Parsed() {
		if *gffFlag == "" {
			fmt.Println(os.Stderr, "Please supply a gff file as input using -gff option.")
		}
		log.Info(fmt.Sprintf("Converting: %q", *gffFlag))
		coverage.Parse(*gffFlag)

	}
	if mergeCommand.Parsed() {
		if *intervalsFlag == "" {
			fmt.Println(os.Stderr, "Please supply a bed file as input using -intervals option.")
		}
		if *formatFlag != "bed" && *formatFlag != "gff" {
			fmt.Println(os.Stderr, "Please supply a valid format (bed|gff) using -format option.")
			os.Exit(2)
		}
		expand := 20
		if *expandFlag != "" {
			expand, _ = strconv.Atoi(*expandFlag)
		}
		log.Info(fmt.Sprintf("Merging intervals: %q", *intervalsFlag))
		coverage.Merge(*intervalsFlag, expand, *formatFlag)

	}
}

func serve() {
	log.Printf("main: starting HTTP server")
	srv := coverage.Serve()
	//log.Printf("main: serving for 60 seconds")
	if err := srv.ListenAndServe(); err != nil {
		// cannot panic, because this probably is an intentional close
		log.Warn(fmt.Sprintf("Httpserver: ListenAndServe() error: %s", err))
	}

	//log.Printf("main: stopping HTTP server")
	//
	//if err := srv.Shutdown(nil); err != nil {
	//	panic(err) // failure/timeout shutting down the server gracefully
	//}
	//
	//log.Printf("main: done. exiting")
}
