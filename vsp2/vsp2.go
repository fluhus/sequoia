// Extracts enriched species from the VSP2 table.
package main

import (
	"fmt"
	"lab/common"
	"regexp"

	"github.com/fluhus/gostuff/csvdec"
	"github.com/fluhus/gostuff/jio"
)

const (
	inFile  = "vsp2_species.csv"
	outFile = "vsp2_species.json"
	// inFile  = "vsp2_species_official.csv"
	// outFile = "vsp2_species2.json"
)

func main() {
	m := map[string][]string{}
	noIDs := 0
	for e, err := range csvdec.FileHeader[entry](inFile, nil) {
		common.Die(err)
		if len(e.IDs) == 0 {
			noIDs++
			// continue
		}
		m[e.Name] = e.IDs
	}
	fmt.Println(len(m), "no IDs:", noIDs)
	if true {
		jio.Write(outFile, m)
	}
}

type entry struct {
	Name string   `csvdec:"Reporting Name"`
	IDs  []string `csvdec:"Genome Length Reference,ParseIDs"`
}

var accsRE = regexp.MustCompile(`[;,]\s*`)

func (e entry) ParseIDs(s string) ([]string, error) {
	if s == "" || s == "-" {
		return nil, nil
	}
	return accsRE.Split(s, -1), nil
}
