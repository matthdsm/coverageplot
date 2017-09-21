/*
The coverageplot command performs a coverage analysis on sequencing data.
You can use the d3.js based interface to visualise the results of the analysis
in an exon centric way.
*/
package main

import (
	"flag"
	"fmt"
	coverage "gitlab.cmgg.be/cmgg/coverageplot"
	"log"
	"os"
	"strconv"
	"time"
)

func main() {
	analyseCommand := flag.NewFlagSet("analyse", flag.ExitOnError)
	analyseIntervalsFlag := analyseCommand.String("intervals", "", "Input the sorted and merged interval file.")
	bamFlag := analyseCommand.String("bam", "", "Input bam file, a sorted and indexed file is required.")

	serveCommand := flag.NewFlagSet("serve", flag.ExitOnError)

	intervalsCommand := flag.NewFlagSet("intervals", flag.ExitOnError)
	gffFlag := intervalsCommand.String("gff", "", "Input gff file.")

	mergeCommand := flag.NewFlagSet("merge", flag.ExitOnError)
	intervalsFlag := mergeCommand.String("intervals", "", "Input the sorted interval file.")
	expandFlag := mergeCommand.String("expand", "", "Expand the interval with # bp on each side")

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
	default:
		fmt.Fprintf(os.Stderr, "%q is not valid command.\n", os.Args[1])
		os.Exit(2)
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
		fmt.Fprintf(os.Stderr, "Analysing bam file %q using intervals %q\n", *bamFlag, *analyseIntervalsFlag)
		coverage.Analyse(*bamFlag, *analyseIntervalsFlag)
	}

	if serveCommand.Parsed() {
		serve()
	}

	if intervalsCommand.Parsed() {
		if *gffFlag == "" {
			fmt.Println(os.Stderr, "Please supply a gff file as input using -gff option.")
		}
		fmt.Fprintf(os.Stderr, "Converting: %q\n", *gffFlag)
		coverage.Parse(*gffFlag)

	}
	if mergeCommand.Parsed() {
		if *intervalsFlag == "" {
			fmt.Println(os.Stderr, "Please supply a bed file as input using -intervals option.")
		}
		expand := 20
		if *expandFlag != "" {
			expand, _ = strconv.Atoi(*expandFlag)
		}
		fmt.Fprintf(os.Stderr, "Merging intervals: %q\n", *intervalsFlag)
		coverage.Merge(*intervalsFlag, expand)

	}
}

func serve() {
	log.Printf("main: starting HTTP server")
	srv := coverage.Serve()
	log.Printf("main: serving for 60 seconds")
	time.Sleep(60 * time.Second)

	log.Printf("main: stopping HTTP server")

	if err := srv.Shutdown(nil); err != nil {
		panic(err) // failure/timeout shutting down the server gracefully
	}

	log.Printf("main: done. exiting")
}
