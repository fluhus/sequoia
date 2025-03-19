// Extracts important data from the raw refseq dump.
package main

import (
	"bytes"
	"fmt"
	"lab/common"
	"lab/config"
	"lab/genbankplus"
	"runtime/debug"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/biostuff/formats/genbank"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/ptimer"
	"github.com/fluhus/gostuff/sets"
)

const (
	mustHaveSequence = false

	dir            = config.DataDir + "/refseq/"
	pgrFile        = dir + "viral.1.genomic.gbff.gz"
	missingPGRFile = dir + "vsp2_missing.pgr"
	outMeta        = dir + "viral.1.genomic.json"
	outSeqs        = dir + "viral.1.genomic.fasta.zst"
	outTID2NC      = dir + "tid2nc.json"
	outNC2TID      = dir + "nc2tid.json"
)

func main() {
	debug.SetGCPercent(20)

	fmt.Println("Reading data")
	fseq, err := aio.Create(outSeqs)
	common.Die(err)
	var ee []genbankplus.GenBankPlus
	pt := ptimer.New()
	noSeq := 0
	foundAccs := sets.Set[string]{}
	foundTIDs := sets.Set[string]{}
	for _, file := range []string{pgrFile, missingPGRFile} {
		for e, err := range genbankplus.Iter(genbank.File(file)) {
			common.Die(err)
			if len(e.Origin) == 0 {
				if mustHaveSequence {
					common.Die(fmt.Errorf("no sequence:\n%v", *e))
				}
				noSeq++
				continue
			}
			if foundAccs.Has(e.Accessions[0]) {
				common.Die(fmt.Errorf("duplicate acc: %q", e.Accessions[0]))
			}
			foundAccs.Add(e.Accessions[0])
			foundTIDs.Add(e.TaxID)
			fa := &fasta.Fasta{
				Name:     []byte(e.Accessions[0]),
				Sequence: bytes.ToUpper([]byte(e.Origin)),
			}
			common.Die(fa.Write(fseq))
			e.Origin = ""
			e.References = nil
			e.Features = nil
			e.DBLink = nil
			ee = append(ee, *e)
			pt.Inc()
		}
	}
	common.Die(fseq.Close())
	pt.Done()
	fmt.Println("Unique taxids:", len(foundTIDs))
	fmt.Println("noSeq", noSeq)

	fmt.Println("Writing metadata")
	pt = ptimer.New()
	common.Die(jio.Write(outMeta, ee))
	pt.Done()

	fmt.Println("Creating mapping")
	pt = ptimer.New()
	tid2nc := map[string][]string{}
	nc2tid := map[string]string{}
	muktzeh := sets.Set[string]{}
	for _, e := range ee {
		// Using the first acc as sequence identifier.
		// tid2nc[e.TaxID] = append(tid2nc[e.TaxID], e.Accs[0])
		tid2nc[e.TaxID] = append(tid2nc[e.TaxID], e.Accessions...)
		for _, acc := range e.Accessions {
			if nc2tid[acc] != "" {
				// common.Die(fmt.Errorf("duplicate acc: %q, tids %q %q",
				// 	acc, nc2tid[acc], e.TaxID))
				muktzeh.Add(acc)
			}
			nc2tid[acc] = e.TaxID
		}
	}
	fmt.Println(len(muktzeh), "accs with multiple tids")
	for acc := range muktzeh {
		delete(nc2tid, acc)
	}
	for k, v := range extraMapping {
		tid2nc[k] = append(tid2nc[k], v)
	}

	fmt.Println(len(tid2nc), "tids")
	fmt.Println(len(nc2tid), "accs")
	common.Die(jio.Write(outTID2NC, tid2nc))
	common.Die(jio.Write(outNC2TID, nc2tid))
	pt.Done()
}

// Manually located stuff.
var extraMapping = map[string]string{
	"2844245": "NC_049948",
	"1986019": "NC_022266",
}
