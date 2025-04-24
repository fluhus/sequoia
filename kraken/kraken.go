// Package kraken parses kraken and bracken reports.
package kraken

import (
	"encoding/csv"
	"fmt"
	"iter"
	"slices"
	"strconv"
	"strings"

	"github.com/fluhus/gostuff/csvdec"
)

func File(file string) iter.Seq2[Entry, error] {
	return func(yield func(Entry, error) bool) {
		hierarchy := make([]string, 0, 10)
		for e, err := range csvdec.File[Entry](file, toTSV) {
			if err != nil {
				yield(e, err)
				return
			}
			level, err := nameToLevel(e.Name)
			if err != nil {
				yield(e, err)
				return
			}
			e.Level = level
			hierarchy = append(hierarchy[:level], strings.Trim(e.Name, " "))
			e.Hierarchy = slices.Clone(hierarchy)
			if !yield(e, nil) {
				return
			}
		}
	}
}

type Entry struct {
	Perc        float64 `csvdec:",ParsePerc"`
	ReadsClade  int
	ReadsDirect int
	Rank        string
	TaxID       int
	Name        string
	Level       int      `csvdec:"-"`
	Hierarchy   []string `csvdec:"-"`
}

func (k Entry) ParsePerc(s string) (float64, error) {
	return strconv.ParseFloat(strings.Trim(s, " "), 64)
}

func toTSV(r *csv.Reader) {
	r.Comma = '\t'
	r.FieldsPerRecord = -1
}

func nameToLevel(name string) (int, error) {
	for i, x := range name {
		if x != ' ' { // Hit a non-space.
			if i%2 != 0 { // Odd number of spaces, not good!
				return 0, fmt.Errorf("bad number of spaces in name: %v", i)
			}
			return i / 2, nil
		}
	}
	return 0, fmt.Errorf("found only spaces in name: %q", name)
}
