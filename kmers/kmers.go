// Extracts kmers from samples.
package main

import (
	"flag"
	"fmt"
	"iter"
	"lab/common"
	"path/filepath"
	"runtime/debug"
	"slices"

	"github.com/fluhus/biostuff/formats/fastq"
	"github.com/fluhus/biostuff/sequtil"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/bits"
	"github.com/fluhus/gostuff/bnry"
	"github.com/fluhus/gostuff/hashx"
	"github.com/fluhus/gostuff/iterx"
	"github.com/fluhus/gostuff/ptimer"
)

const (
	nbits    = 32
	nbytes   = 1 << (nbits - 3)
	mask     = 1<<nbits - 1
	kmerLen  = 31
	maxReads = 10000000
	ratio    = 32
)

var (
	inGlob  = flag.String("i", "", "Input file glob")
	outFile = flag.String("o", "", "Output file")
)

func main() {
	debug.SetGCPercent(20)

	flag.Parse()
	inFiles, err := filepath.Glob(*inGlob)
	common.Die(err)

	found := make([]byte, nbytes)
	found2 := make([]byte, nbytes)
	fmt.Println("nbytes:", nbytes)

	for _, inFile := range inFiles {
		fmt.Println(inFile)
		pt := ptimer.New()
		for fq, err := range iterx.Limit2(fastq.File(inFile), maxReads) {
			common.Die(err)
			for h := range seqKmersHashes(fq.Sequence, kmerLen, ratio) {
				if bits.Get(found, h&mask) == 1 {
					bits.Set1(found2, h&mask)
				} else {
					bits.Set1(found, h&mask)
				}
			}
			pt.Inc()
		}
		pt.Done()
		fmt.Println(bits.Sum(found), bits.Sum(found2))
	}

	fmt.Println("Writing")
	fout, err := aio.Create(*outFile)
	common.Die(err)
	// common.Die(bnry.Write(fout, found))
	ones := slices.Collect(bits.Ones(found2))
	common.ToDiffs(ones)
	common.Die(bnry.Write(fout, ones))
	fout.Close()
}

var h1 = hashx.NewSeed(0)
var h2 = hashx.NewSeed(1)

func seqKmersHashes(seq []byte, k int, p uint64) iter.Seq[uint64] {
	return func(yield func(uint64) bool) {
		for kmer := range sequtil.CanonicalSubsequences(seq, k) {
			if p > 1 && h1.Bytes(kmer)%p != 0 {
				continue
			}
			if !yield(h2.Bytes(kmer)) {
				return
			}
		}
	}
}
