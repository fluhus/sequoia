// Filters output reads from kraken, keeping only ones that mapped
// to specific tax IDs.
package main

import (
	"fmt"
	"lab/common"
	"os"
	"regexp"

	"github.com/fluhus/biostuff/formats/fastq"
	"github.com/fluhus/gostuff/sets"
)

func main() {
	// All "Norwalk virus" descendents.
	taxids := sets.Of("11983",
		"95340", "122928", "122929", "262897", "340017", "1246677",
		"235544", "490039", "490043", "552592", "1160947",
		"1529909", "1529918", "1529924")
	fmt.Fprintln(os.Stderr, len(taxids), "taxids")

	re := regexp.MustCompile(`kraken:taxid\|(\d+)`)
	for fq, err := range fastq.Reader(os.Stdin) {
		common.Die(err)
		m := re.FindSubmatch(fq.Name)
		if m == nil {
			common.Die(fmt.Errorf("no kraken taxid: %q", fq.Name))
		}
		if taxids.Has(string(m[1])) {
			fq.Write(os.Stdout)
		}
	}
}
