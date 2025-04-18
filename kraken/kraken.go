// Converts kraken reports to JSON.
package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"lab/common"
	"strconv"
	"strings"

	"github.com/fluhus/gostuff/csvdec"
	"github.com/fluhus/gostuff/flagx"
	"github.com/fluhus/gostuff/gnum"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/snm"
	"golang.org/x/exp/maps"
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

	hierarchy := make([]string, 0, 10)
	iabnd := map[string]int{}

	for line, err := range csvdec.File[krakenLine](args[0], toTSV) {
		common.Die(err)

		level, err := nameToLevel(line.Name)
		common.Die(err)
		name := strings.Trim(line.Name, " ")
		hierarchy = append(hierarchy[:level], name) // TODO(amit): Verify we didn't skip a level.

		// Change name after updating hierarchy, to keep hierarchy
		// human-readable.
		if *useTaxIDs {
			name = fmt.Sprint(line.TaxID)
		}

		if line.Rank == "U" && *withOthers {
			iabnd["Unclassified"] += line.ReadsDirect
			continue
		}

		inDomain := level == 2
		var inDomainOfInterest bool
		rankOfInterest := "S"
		switch *mode {
		case "vir":
			inDomainOfInterest = len(hierarchy) >= 2 && hierarchy[1] == "Viruses"
		case "viror":
			inDomainOfInterest = len(hierarchy) >= 2 && hierarchy[1] == "Viruses"
		case "virp":
			inDomainOfInterest = len(hierarchy) >= 2 && hierarchy[1] == "Viruses"
		case "virg":
			inDomainOfInterest = len(hierarchy) >= 2 && hierarchy[1] == "Viruses"
			rankOfInterest = "G"
		case "bact":
			inDomainOfInterest = len(hierarchy) >= 3 && hierarchy[2] == "Bacteria"
		case "allsp":
			inDomainOfInterest = true // Report all species.
		case "allgen":
			rankOfInterest = "G"
			inDomainOfInterest = true // Report all genuses.
		}

		if !inDomainOfInterest {
			if *mode == "viror" {
				iabnd["Other"] += line.ReadsDirect
			} else if *withOthers && inDomain {
				iabnd[name] += line.ReadsClade
			}
			continue
		}

		reads := 0
		if rankOfInterest == "S" {
			// If want "S", acccept "S1", "S2"...
			if !strings.HasPrefix(line.Rank, rankOfInterest) {
				continue
			}
			// Report only leaves.
			// if line.ReadsClade != line.ReadsDirect {
			// 	continue
			// }
			reads = line.ReadsDirect
			if reads == 0 {
				continue
			}
		} else {
			if line.Rank != rankOfInterest {
				continue
			}
			reads = line.ReadsClade
		}

		if *mode == "virp" { // Attach phylum to the name.
			if len(hierarchy) < 5 {
				name = fmt.Sprint(name, ",unknown")
			} else {
				name = fmt.Sprint(name, ",", hierarchy[4])
			}
		}
		iabnd[name] += reads
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

type krakenLine struct {
	Perc        float64 `csvdec:",ParsePerc"`
	ReadsClade  int
	ReadsDirect int
	Rank        string
	TaxID       int
	Name        string
}

func (k krakenLine) ParsePerc(s string) (float64, error) {
	return strconv.ParseFloat(strings.Trim(s, " "), 64)
}

func toTSV(r *csv.Reader) {
	r.Comma = '\t'
	r.FieldsPerRecord = -1
}

func nameToLevel(name string) (int, error) {
	for i, x := range name {
		if x != ' ' { // Hit a non-space.
			if i%2 != 0 { // Odd number of spaces, not good!
				return 0, fmt.Errorf("bad number of spaces in name: %v", i)
			}
			return i / 2, nil
		}
	}
	return 0, fmt.Errorf("found only spaces in name: %q", name)
}
