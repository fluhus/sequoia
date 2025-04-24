// Converts kraken reports to JSON.
package main

import (
	"flag"
	"fmt"
	"lab/common"
	"lab/kraken"
	"strings"

	"github.com/fluhus/gostuff/flagx"
	"github.com/fluhus/gostuff/gnum"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/snm"
	"golang.org/x/exp/maps"
)

const (
	verbose = false
)

var (
	useTaxIDs     = flag.Bool("tid", false, "Use taxon IDs rather than taxon names")
	withOthers    = flag.Bool("others", false, "Include other domains")
	keepReadCount = flag.Bool("readcount", false, "Output read counts without normalizing to 1")
	mode          = flagx.OneOf("m", "vir", "Reporting mode",
		"vir", "virp", "virg", "bact", "viror", "allsp", "allgen")
)

func main() {
	flag.Parse()
	args := flag.Args()

	if verbose {
		fmt.Println("Running on file:", args[0])
		fmt.Println("Mode:", *mode)
	}

	iabnd := map[string]int{}
	i := 0

	for e, err := range kraken.File(args[0]) {
		common.Die(err)
		i++

		name := strings.Trim(e.Name, " ")
		if *useTaxIDs {
			name = fmt.Sprint(e.TaxID)
		}

		if e.Rank == "U" && *withOthers {
			iabnd["Unclassified"] += e.ReadsDirect
			continue
		}

		inDomain := e.Level == 2
		var inDomainOfInterest bool
		rankOfInterest := "S"
		if verbose {
			fmt.Println(e.Hierarchy)
		}
		switch *mode {
		case "vir":
			inDomainOfInterest = len(e.Hierarchy) >= 2 && e.Hierarchy[1] == "Viruses"
		case "viror":
			inDomainOfInterest = len(e.Hierarchy) >= 2 && e.Hierarchy[1] == "Viruses"
		case "virp":
			inDomainOfInterest = len(e.Hierarchy) >= 2 && e.Hierarchy[1] == "Viruses"
		case "virg":
			inDomainOfInterest = len(e.Hierarchy) >= 2 && e.Hierarchy[1] == "Viruses"
			rankOfInterest = "G"
		case "bact":
			inDomainOfInterest = len(e.Hierarchy) >= 3 && e.Hierarchy[2] == "Bacteria"
		case "allsp":
			inDomainOfInterest = true // Report all species.
		case "allgen":
			rankOfInterest = "G"
			inDomainOfInterest = true // Report all genuses.
		}

		if !inDomainOfInterest {
			if *mode == "viror" {
				iabnd["Other"] += e.ReadsDirect
			} else if *withOthers && inDomain {
				iabnd[name] += e.ReadsClade
			}
			continue
		}

		reads := 0
		if rankOfInterest == "S" {
			// If want "S", acccept "S1", "S2"...
			if !strings.HasPrefix(e.Rank, rankOfInterest) {
				continue
			}
			// Report only leaves.
			// if line.ReadsClade != line.ReadsDirect {
			// 	continue
			// }
			reads = e.ReadsDirect
			if reads == 0 {
				continue
			}
		} else {
			if e.Rank != rankOfInterest {
				continue
			}
			reads = e.ReadsClade
		}

		if *mode == "virp" { // Attach phylum to the name.
			if len(e.Hierarchy) < 5 {
				name = fmt.Sprint(name, ",unknown")
			} else {
				name = fmt.Sprint(name, ",", e.Hierarchy[4])
			}
		}
		iabnd[name] += reads
	}
	if verbose {
		fmt.Println("Parsed", i, "lines")
	}

	sum := gnum.Sum(maps.Values(iabnd))
	if *keepReadCount {
		sum = 1
	}
	abnd := snm.MapToMap(iabnd, func(k string, v int) (string, float64) {
		return k, float64(v) / float64(sum)
	})
	common.Die(jio.Write(args[1], abnd))
}
