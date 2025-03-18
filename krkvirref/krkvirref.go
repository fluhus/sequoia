// Extract sequences from refseq for kraken's taxa.
package main

import (
	"fmt"
	"lab/common"
	"runtime/debug"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/sets"
)

const (
	krkFile    = "krk_viral.json"
	refseqFile = "/dfs7/whitesonlab/alavon/Data/refseq/viral.1.genomic.fasta.zst"
	outFile    = "/dfs7/whitesonlab/alavon/Data/refseq/kraken.fa.zst"
)

func main() {
	debug.SetGCPercent(20)

	accs := sets.Set[string]{}
	{
		m := map[string][]string{}
		common.Die(jio.Read(krkFile, &m))
		for _, v := range m {
			accs.Add(v...)
		}
	}
	fmt.Println(len(accs))

	f, err := aio.Create(outFile)
	common.Die(err)
	defer f.Close()
	found := 0
	for fa, err := range fasta.File(refseqFile) {
		common.Die(err)
		if !accs.Has(string(fa.Name)) {
			continue
		}
		common.Die(fa.Write(f))
		found++
	}
	fmt.Println("Found", found)
}
