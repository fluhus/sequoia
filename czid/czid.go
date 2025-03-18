// Parses CZID output files.
package main

import (
	"fmt"
	"lab/common"
	"os"
	"strconv"
	"strings"

	"github.com/fluhus/gostuff/csvdec"
	"github.com/fluhus/gostuff/gnum"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/sets"
	"golang.org/x/exp/maps"
)

func main() {
	if len(os.Args) != 3 {
		fmt.Println("Usage: czid IN_FILE OUT_FILE")
		os.Exit(1)
	}
	abnd := map[string]float64{}
	bactNames := sets.Set[string]{}
	for e, err := range csvdec.FileHeader[Entry](os.Args[1], nil) {
		common.Die(err)
		if e.NtRPM == 0 {
			continue
		}
		if e.TaxLevel != 2 {
			continue
		}
		if e.Category == "bacteria" {
			bactNames.Add(e.Name)
		}
		if e.Category != "viruses" {
			continue
		}
		abnd[e.Name] += e.NtRPM
	}
	fmt.Println("Bacterial species:", len(bactNames))
	fmt.Println("Viral species:    ", len(abnd))

	sum := gnum.Sum(maps.Values(abnd))
	for k := range abnd {
		abnd[k] /= sum
	}
	jio.Write(os.Args[2], abnd)
}

func (Entry) ParseTaxIDs(s string) ([]int, error) {
	if s == "" || s[0] != '[' || s[len(s)-1] != ']' {
		return nil, fmt.Errorf("could not parse tax IDs: %q", s)
	}
	parts := strings.Split(s[1:len(s)-1], ", ")
	result := make([]int, 0, len(parts))
	for _, p := range parts {
		i, err := strconv.Atoi(p)
		if err != nil {
			return nil, err
		}
		result = append(result, i)
	}
	return result, nil
}

type Entry struct {
	Name          string
	TaxLevel      int     `csvdec:"tax_level"`
	NtRPM         float64 `csvdec:"nt_rpm,allowempty"`
	NrRPM         float64 `csvdec:"nr_rpm,allowempty"`
	SpeciesTaxIDs []int   `csvdec:"species_tax_ids,allowempty,ParseTaxIDs"`
	Category      string
}
