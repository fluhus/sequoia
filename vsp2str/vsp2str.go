// String matching between vsp2 species names and kraken results.
package main

import (
	"fmt"
	"lab/common"
	"math"
	"os"
	"strings"

	"github.com/fluhus/biostuff/align"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/sets"
	"github.com/fluhus/gostuff/snm"
	"golang.org/x/exp/maps"
)

func main() {
	vsp2, err := readVSP2Names()
	common.Die(err)
	fmt.Println(len(vsp2), "vsp2 species")
	// fmt.Println(vsp2[:5])

	idx := map[string][]int{}
	for i, name := range vsp2 {
		for ss := range substrs(name, 4) {
			idx[ss] = append(idx[ss], i)
		}
	}
	fmt.Println(len(idx), "substrs in index")

	idx = filterIndex(idx, isqrt(len(vsp2)))
	fmt.Println(len(idx), "substrs in index (post filtering)")

	queries, err := readAbundanceNames(os.Args[1])
	common.Die(err)
	fmt.Println(len(queries), "queries")
	for _, q := range queries {
		found := searchIndex(idx, q)
		names := snm.SliceToSlice(found, func(i int) string { return vsp2[i] })
		names = snm.FilterSlice(names, func(s string) bool {
			ln := max(len(q), len(s))
			_, x := align.Global([]byte(q), []byte(s), align.Levenshtein)
			return x > -float64(ln)/5
		})
		if len(names) == 0 {
			continue
		}
		fmt.Println(q, names)
	}
}

func readVSP2Names() ([]string, error) {
	m := map[string][]string{}
	if err := jio.Read("vsp2_species.json", &m); err != nil {
		return nil, err
	}
	keys := maps.Keys(m)
	for i := range keys {
		keys[i] = strings.ToLower(keys[i])
	}
	return keys, nil
}

func readAbundanceNames(file string) ([]string, error) {
	m := map[string]float64{}
	if err := jio.Read(file, &m); err != nil {
		return nil, err
	}
	keys := maps.Keys(m)
	for i := range keys {
		keys[i] = strings.ToLower(keys[i])
	}
	return keys, nil
}

func substrs(s string, n int) sets.Set[string] {
	ss := sets.Set[string]{}
	for i := range len(s) - n + 1 {
		ss.Add(s[i : i+n])
	}
	return ss
}

func isqrt(i int) int {
	return int(math.Round(math.Sqrt(float64(i))))
}

func filterIndex(idx map[string][]int, mx int) map[string][]int {
	idx2 := map[string][]int{}
	for k, v := range idx {
		if len(v) > mx {
			// fmt.Println("-", k)
			continue
		}
		idx2[k] = v
	}
	return idx2
}

func searchIndex(idx map[string][]int, s string) []int {
	found := sets.Set[int]{}
	for ss := range substrs(s, 4) {
		for _, i := range idx[ss] {
			found.Add(i)
		}
	}
	return maps.Keys(found)
}
