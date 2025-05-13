// Sums up coverage histograms and groups them by sample type.
package main

import (
	"fmt"
	"iter"
	"lab/common"
	"lab/luna"
	"math"
	"os"
	"path/filepath"
	"regexp"
	"runtime/debug"
	"slices"
	"strings"

	"github.com/fluhus/gostuff/gnum"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/snm"
	"golang.org/x/exp/maps"
)

const (
	meanType = linMean
	batch    = 1
)

func main() {
	debug.SetGCPercent(20)

	inDir := os.Args[1]
	files, err := sortFiles(filepath.Join(inDir, "*.cov.json"))
	common.Die(err)
	fmt.Println(maps.Keys(files))
	fmt.Println("Mean type:", meanType)

	for spc, it := range iterFiles(files) {
		fmt.Println(spc)
		sumByGroup := map[string][]int{}
		if batch == 1 {
			for g := range luna.GroupNameMapping {
				sumByGroup[g] = nil
			}
		}
		logByGroup := map[string][]float64{}
		nz := map[string]float64{}
		nSamples := 0
		for sm, err := range it {
			common.Die(err)
			if batch == 1 && strings.HasPrefix(sm.sample, "Undetermined") {
				continue
			}
			if batch == 2 && !batch2RE.MatchString(sm.sample) {
				continue
			}
			nSamples++
			fmt.Println("--", sm.sample, len(sm.cov))
			if batch == 1 {
				sm.sample = luna.FixName(sm.sample)
			}
			g := sampleGroup(sm.sample)
			if meanType == geoMean {
				logByGroup[g] = addLog(logByGroup[g], sm.cov)
			} else {
				sumByGroup[g] = add(sumByGroup[g], sm.cov)
			}
			nz[sm.sample] = nonZeroRatio(sm.cov)
		}
		fmt.Printf("^^ (%v samples)\n", nSamples)

		// Divide by number of samples for mean.
		switch meanType {
		case linMean:
			for _, v := range sumByGroup {
				for i, x := range v {
					v[i] = gnum.Idiv(x, nSamples)
				}
			}
		case geoMean:
			for k, v := range logByGroup {
				sumByGroup[k] = snm.SliceToSlice(v, func(f float64) int {
					// -1 because we added 1 before log.
					return int(math.Round(math.Exp(f/float64(nSamples)))) - 1
				})
			}
		}

		if batch == 1 && len(sumByGroup) != 4 { // Assertion
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
	if len(files) == 0 {
		return nil, fmt.Errorf("no files found")
	}
	re := regexp.MustCompile(`([^/]+)__(.+)\.cov\.json`)
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
	dst = grow(dst, len(src))
	for i, x := range src {
		dst[i] += x
	}
	return dst
}

func addLog(dst []float64, src []int) []float64 {
	dst = grow(dst, len(src))
	for i, x := range src {
		dst[i] += math.Log(float64(x + 1))
	}
	return dst
}

func grow[T any](a []T, n int) []T {
	if len(a) >= n {
		return a
	}
	a = slices.Grow(a, n-len(a))
	var zero T
	for len(a) < n {
		a = append(a, zero)
	}
	return a
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
		m[k] = grow(v, mx)
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

var batch2RE = regexp.MustCompile(`_([^_]+)_(SOL|INF)_`)

func sampleGroup(s string) string {
	if batch == 1 {
		return luna.SampleGroup(s)
	}
	if batch == 2 {
		m := batch2RE.FindStringSubmatch(s)
		if m == nil {
			panic(fmt.Sprintf("bad sample name: %q", s))
		}
		return m[1] + "_" + m[2]
	}
	panic(fmt.Sprintf("unreachable, batch=%v", batch))
}

// Gurads batch to be 1 or 2.
func _() {
	var x [2]struct{}
	_ = x[batch-1]
}

const (
	noMean = iota
	linMean
	geoMean
)
