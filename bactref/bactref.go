// Creates fasta files from a genbank file.
package main

import (
	"fmt"
	"lab/common"
	"lab/genbankplus"
	"os"
	"path/filepath"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/biostuff/formats/genbank"
	"github.com/fluhus/gostuff/aio"
)

func main() {
	if len(os.Args) != 3 {
		fmt.Printf("Usage: %v IN_PGR OUT_DIR\n", filepath.Base(os.Args[0]))
		os.Exit(1)
	}
	inPGR, outDir := os.Args[1], os.Args[2]

	common.Die(os.MkdirAll(outDir, 0o755))

	for e, err := range genbankplus.Iter(genbank.File(inPGR)) {
		common.Die(err)
		fmt.Println(e.TaxID, e.Accessions, len(e.Origin))
		acc := e.Accessions[0]
		fa := &fasta.Fasta{Name: []byte(acc), Sequence: []byte(e.Origin)}
		f, err := aio.Create(filepath.Join(outDir, acc+".fasta"))
		common.Die(err)
		common.Die(fa.Write(f))
		f.Close()
	}
}
