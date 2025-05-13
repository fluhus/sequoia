// Creates rarefaction curves.
package main

import (
	"encoding/csv"
	"fmt"
	"iter"
	"lab/common"
	"os"
	"path/filepath"
	"regexp"
	"runtime/debug"
	"strconv"
	"strings"

	"github.com/fluhus/biostuff/rarefy"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/csvdec"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/ptimer"
	"github.com/fluhus/gostuff/sets"
	"golang.org/x/exp/maps"
)

const (
	useAllSpecies = false
	nperms        = 10
	rstep         = 1000
)

func main() {
	debug.SetGCPercent(20)

	// ezpprof.Start("rarefy/default.pgo")
	// ezpprof.Start("profile_data")
	// defer ezpprof.Stop()

	files, _ := filepath.Glob(os.Args[1])
	if len(files) == 0 {
		common.Die(fmt.Errorf("no files found"))
	}
	// files = files[:4]
	fmt.Println("Found", len(files), "files")

	fmt.Println("Loading abundances")
	var names []string
	var counts []map[string]int
	pt := ptimer.New()
	for _, f := range files {
		m, err := readFile(f)
		common.Die(err)
		counts = append(counts, m)
		name := strings.TrimSuffix(filepath.Base(f), ".krk.txt")
		names = append(names, name)
		pt.Inc()
	}
	pt.Done()

	fmt.Println("Saving")
	common.Die(writeCSV(toTable(names, counts), os.Args[2]))

	fmt.Println("Rarefying")
	pt = ptimer.New()
	rrf := map[string][][]int{}
	for i := range names {
		xx, yy := rarefy.Rarefy(maps.Values(counts[i]), rstep, nperms)
		rrf[names[i]] = [][]int{xx, yy}
		pt.Inc()
	}
	pt.Done()
	common.Die(jio.Write("rrf.json", rrf))
}

func readFile(file string) (map[string]int, error) {
	inViruses := false
	levelRE := regexp.MustCompile(`^ ? ?\S`)
	domainRE := regexp.MustCompile(`^    \S`)
	m := map[string]int{}

	for line, err := range csvdec.File[krakenLine](file, toTSV) {
		if err != nil {
			return nil, err
		}
		if line.Rank == "U" {
			m["Other"] += line.ReadsDirect
			continue
		}
		if levelRE.MatchString(line.Name) {
			inViruses = line.Name == "  Viruses"
		}
		if !useAllSpecies && !inViruses {
			if domainRE.MatchString(line.Name) && line.ReadsClade > 0 {
				// name := strings.Trim(line.Name, " ")
				m["Other"] += line.ReadsClade
			}
			continue
		}
		if line.Rank != "S" {
			continue
		}
		name := strings.Trim(line.Name, " ")
		m[name] += line.ReadsClade
	}
	return m, nil
}

type krakenLine struct {
	Perc        float64 `csvdec:",ParsePerc"`
	ReadsClade  int
	ReadsDirect int
	Rank        string
	TaxID       int
	Name        string
}

func (k krakenLine) ParsePerc(s string) (float64, error) {
	return strconv.ParseFloat(strings.Trim(s, " "), 64)
}

func toTSV(r *csv.Reader) {
	r.Comma = '\t'
	r.FieldsPerRecord = -1
}

func toTable(names []string, vals []map[string]int) iter.Seq[[]string] {
	return func(yield func([]string) bool) {
		colSet := sets.Set[string]{}
		for _, v := range vals {
			sets.AddKeys(colSet, v)
		}
		cols := maps.Keys(colSet)

		var row []string
		row = append(row, "name")
		row = append(row, cols...)
		if !yield(row) {
			return
		}
		for i, n := range names {
			row = row[:0]
			row = append(row, n)
			for _, c := range cols {
				row = append(row, fmt.Sprint(vals[i][c]))
				// if vals[i][c] == 0 {
				// 	panic(fmt.Sprintf("waaat %q %q", n, c))
				// }
			}
			if !yield(row) {
				return
			}
		}
	}
}

func writeCSV(rows iter.Seq[[]string], file string) error {
	f, err := aio.Create(file)
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	defer w.Flush()
	for row := range rows {
		if err := w.Write(row); err != nil {
			return err
		}
	}
	return nil
}
