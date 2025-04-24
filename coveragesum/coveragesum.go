// Sums up coverage histograms and groups them by sample type.
package main

import (
	"fmt"
	"iter"
	"lab/common"
	"lab/luna"
	"os"
	"path/filepath"
	"regexp"
	"runtime/debug"
	"strings"

	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/snm"
	"golang.org/x/exp/maps"
)

func main() {
	debug.SetGCPercent(20)

	inDir := os.Args[1]
	files, err := sortFiles(filepath.Join(inDir, "*.cov.json"))
	common.Die(err)
	fmt.Println(maps.Keys(files))

	sampleToName := sampleNumToName()

	for spc, it := range iterFiles(files) {
		fmt.Println(spc)
		sumByGroup := map[string][]int{}
		for g := range luna.GroupNameMapping {
			sumByGroup[g] = nil
		}
		nz := map[string]float64{}
		for sm, err := range it {
			common.Die(err)
			fmt.Println("--", sm.sample, len(sm.cov))
			name := sampleToName[sm.sample]
			if strings.HasPrefix(name, "Undetermined") {
				continue
			}
			name = luna.FixName(sampleToName[sm.sample])
			// if name == "" {
			// 	common.Die(fmt.Errorf("bad name: %q", sm.sample))
			// }
			// g := luna.SampleGroup(luna.FixName(name))
			g := luna.SampleGroup(name)
			sumByGroup[g] = add(sumByGroup[g], sm.cov)
			nz[name] = nonZeroRatio(sm.cov)
		}
		if len(sumByGroup) != 4 { // Assertion
			panic(fmt.Sprintf("bad length: %v", len(sumByGroup)))
		}
		equalizeLens(sumByGroup)
		common.Die(jio.Write(filepath.Join(inDir, spc+".covs.json"), sumByGroup))
		common.Die(jio.Write(filepath.Join(inDir, spc+".nz.json"), nz))
	}
}

func sortFiles(glob string) (map[string]map[string]string, error) {
	files, err := filepath.Glob(glob)
	if err != nil {
		return nil, err
	}
	re := regexp.MustCompile(`([^/]+)\.(\d+)\.cov\.json`)
	m := snm.NewDefaultMap(func(s string) map[string]string {
		return map[string]string{}
	})
	for _, f := range files {
		mch := re.FindAllStringSubmatch(f, -1)
		if mch == nil {
			return nil, fmt.Errorf("file name does not match pattern: %q", f)
		}
		species, sample := mch[0][1], mch[0][2]
		m.Get(species)[sample] = f
	}
	return m.M, nil
}

func iterFiles(files map[string]map[string]string) iter.Seq2[string, iter.Seq2[sampleCoverage, error]] {
	return func(yield func(string, iter.Seq2[sampleCoverage, error]) bool) {
		for species, f := range files {
			if !yield(species, iterSpecies(f)) {
				return
			}
		}
	}
}

type sampleCoverage struct {
	sample string
	cov    []int
}

func iterSpecies(files map[string]string) iter.Seq2[sampleCoverage, error] {
	return func(yield func(sampleCoverage, error) bool) {
		for sample, file := range files {
			m := map[string][]int{}
			if err := jio.Read(file, &m); err != nil {
				yield(sampleCoverage{}, err)
				return
			}
			if len(m) == 0 {
				continue
			}
			if len(m) != 1 {
				yield(sampleCoverage{},
					fmt.Errorf("too many species in input json: %v", len(m)))
				return
			}
			for _, data := range m {
				if !yield(sampleCoverage{sample, data}, nil) {
					return
				}
			}
		}
	}
}

func add(dst, src []int) []int {
	if len(dst) < len(src) {
		dst = append(dst, make([]int, len(src)-len(dst))...)
	}
	for i, x := range src {
		dst[i] += x
	}
	return dst
}

func sampleNumToName() map[string]string {
	m := map[string]string{}
	for i, s := range luna.SampleOrdering {
		m[fmt.Sprint(i+1)] = s
	}
	return m
}

func wrapSlices(m map[string][]int) map[string]any {
	mm := map[string]any{}
	for k, v := range m {
		mm[k] = map[string]any{
			"type": "data",
			"data": v,
		}
	}
	return mm
}

func equalizeLens(m map[string][]int) {
	mx := 0
	for _, v := range m {
		mx = max(mx, len(v))
	}
	if mx == 0 {
		return
	}
	for k, v := range m {
		if len(v) == 0 || len(v) == mx {
			continue
		}
		dif := mx - len(v)
		m[k] = append(v, make([]int, dif)...)
	}
}

func nonZeroRatio(a []int) float64 {
	if len(a) == 0 {
		return 0
	}
	nz := 0
	for _, x := range a {
		if x != 0 {
			nz++
		}
	}
	return float64(nz) / float64(len(a))
}
