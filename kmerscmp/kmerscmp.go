package main

import (
	"encoding/json"
	"fmt"
	"lab/common"
	"lab/config"
	"path/filepath"
	"runtime/debug"
	"strings"

	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/bits"
	"github.com/fluhus/gostuff/bnry"
	"github.com/fluhus/gostuff/ptimer"
)

const (
	inDir = config.WSDataDir + "/kmers"

	nbits    = 32
	nbytes   = 1 << (nbits - 3)
	minCount = 6
)

func main() {
	debug.SetGCPercent(20)

	files, err := filepath.Glob(inDir + "/*.kmers")
	common.Die(err)
	// files = files[:3]
	fmt.Println(len(files), "files")

	fmt.Println("Filtering kmers")
	b1 := make([]byte, nbytes)
	// b2 := make([]byte, nbytes)
	cnt := map[int]int{}
	pt := ptimer.NewFunc(func(i int) string {
		return fmt.Sprintf("%v (%v, %v)", i,
			toMillions(bits.Sum(b1)),
			// toMillions(bits.Sum(b2)),
			toMillions(len(cnt)),
		)
	})
	for _, f := range files {
		kmers, err := readFile(f)
		common.Die(err)
		for _, kmer := range kmers {
			if bits.Get(b1, kmer) == 1 {
				// bits.Set(b2, kmer, true)
				cnt[kmer]++
			} else {
				bits.Set1(b1, kmer)
			}
		}
		pt.Inc()
	}
	pt.Done()

	fmt.Println("All kmers:", bits.Sum(b1))

	table := map[int][]int8{}
	for k, v := range cnt {
		if v >= (minCount - 1) {
			table[k] = make([]int8, len(files))
		}
	}
	fmt.Println("WL kmers: ", len(table))

	fmt.Println("Collecting kmers")
	var names []string
	pt = ptimer.New()
	for i, f := range files {
		kmers, err := readFile(f)
		common.Die(err)
		base := strings.TrimSuffix(filepath.Base(f), ".kmers")
		names = append(names, base)
		for _, kmer := range kmers {
			// if !bits.Get(b2, kmer) {
			// 	continue
			// }
			row, ok := table[kmer]
			if !ok {
				continue
			}
			row[i] = 1
		}
		pt.Inc()
	}
	pt.Done()

	fmt.Println("Writing")
	fout, err := aio.Create("kmers.json")
	common.Die(err)
	j := json.NewEncoder(fout)
	common.Die(j.Encode(names))
	for k, v := range table {
		common.Die(j.Encode([]any{k, v}))
	}
	fout.Close()
}

func readFile(file string) ([]int, error) {
	f, err := aio.Open(file)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	var a []int
	err = bnry.Read(f, &a)
	if err != nil {
		return nil, err
	}
	common.FromDiffs(a)
	return a, nil
}

func toMillions(a int) string {
	return fmt.Sprintf("%.1fM", float64(a)/1000000)
}
