// Creates isolated fastas of enriched species genomes.
package main

import (
	"fmt"
	"lab/common"
	"lab/lazy"
	"strings"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/ptimer"
	"github.com/fluhus/gostuff/sets"
	"github.com/fluhus/gostuff/snm"
)

const (
	outFile = "../data/taxref/%.fasta"
	// ncbiFile = "../../../Data/ncbi/sequence.fasta"
	ncbiFile = "/dfs7/whitesonlab/alavon/Data/refseq/viral.1.genomic.fasta.zst"
)

func main() {
	nc2tax := map[string]string{}
	common.Die(jio.Read("taxref.json", &nc2tax))

	w := snm.NewDefaultMap(func(s string) *lazy.Writer {
		return lazy.Create(strings.ReplaceAll(outFile, "%", s), 1<<20)
	})
	found := sets.Set[string]{}
	pt := ptimer.NewFunc(func(i int) string {
		return fmt.Sprintf("%d (%d files)", i, len(w.M))
	})
	for fa, err := range fasta.File(ncbiFile) {
		common.Die(err)
		pt.Inc()
		tax := nc2tax[string(fa.Name)]
		if tax == "" {
			continue
		}
		common.Die(fa.Write(w.Get(string(fa.Name))))
		found.Add(string(fa.Name))
	}
	for _, ww := range w.M {
		common.Die(ww.Close())
	}
	pt.Done()

	for s := range found {
		delete(nc2tax, s)
	}
	fmt.Println("Not found:", nc2tax)
}
