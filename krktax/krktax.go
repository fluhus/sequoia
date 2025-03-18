// Extracts species metadata from kraken's dataset.
package main

import (
	"fmt"
	"iter"
	"lab/common"
	"os"
	"path/filepath"
	"regexp"

	"github.com/fluhus/gostuff/iterx"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/sets"
	"github.com/fluhus/gostuff/snm"
	"golang.org/x/exp/maps"
)

func main() {
	dbbase, outFile := os.Args[1], os.Args[2]
	seqIDFile := filepath.Join(dbbase, "seqid2taxid.map")
	taxFile := filepath.Join(dbbase, "ktaxonomy.tsv")

	tax, err := readTaxFile(taxFile)
	common.Die(err)
	fmt.Println(len(tax), "taxonomy entries")

	badTIDs := 0
	for e, err := range readSeqIDs(seqIDFile) {
		common.Die(err)
		if tax[e.tid] == nil {
			badTIDs++
			continue
			// common.Die(fmt.Errorf("tid not in taxonomy: %q", e.tid))
		}
		tax[e.tid].Accs = append(tax[e.tid].Accs, e.acc)
	}
	fmt.Println(badTIDs, "TIDs without a taxonomy entry")

	for _, t := range tax {
		if len(t.Accs) <= 1 {
			continue
		}
		t.Accs = snm.Sorted(maps.Keys(sets.Of(t.Accs...)))
	}

	lens := snm.Sorted(snm.SliceToSlice(maps.Values(tax), func(t *taxEntry) int {
		return len(t.Accs)
	}))
	fmt.Println(common.NTiles(100, lens))

	common.Die(jio.Write(outFile, tax))
}

type seqIDEntry struct {
	tid string
	acc string
}

func readSeqIDs(file string) iter.Seq2[seqIDEntry, error] {
	return func(yield func(seqIDEntry, error) bool) {
		re := regexp.MustCompile(`^kraken:taxid\|(\d+)\|(\S+?)(?:\.\d+)?\s+(\S+)$`)
		re2 := regexp.MustCompile(`^(\S+)\s+(\S+)$`)
		i := 0
		for line, err := range iterx.LinesFile(file) {
			if err != nil {
				yield(seqIDEntry{}, err)
				return
			}
			i++
			m := re.FindStringSubmatch(line)
			if m != nil {
				if m[1] != m[3] {
					yield(seqIDEntry{},
						fmt.Errorf("line #%d parts 1 and 3 mismatch: %s",
							i, line))
					return
				}
				if !yield(seqIDEntry{tid: m[1], acc: m[2]}, nil) {
					return
				}
				continue
			}
			m = re2.FindStringSubmatch(line)
			if m != nil {
				if false {
					if !yield(seqIDEntry{tid: m[2], acc: m[1]}, nil) {
						return
					}
				}
				continue
			}
			yield(seqIDEntry{},
				fmt.Errorf("line #%d does not match pattern: %s",
					i, line))
			return
		}
	}
}

type taxEntry struct {
	Name      string
	Level     string
	ParentTID string
	ChildTIDs []string `json:",omitempty"`
	Accs      []string `json:",omitempty"`
}

func readTaxFile(file string) (map[string]*taxEntry, error) {
	m := map[string]*taxEntry{}
	splitter := regexp.MustCompile(`\s*\|\s*`)
	for line, err := range iterx.LinesFile(file) {
		if err != nil {
			return nil, err
		}
		parts := splitter.Split(line, -1)
		if len(parts) != 5 {
			return nil, fmt.Errorf("bad number of parts: %v, want 5: %q",
				len(parts), line)
		}
		taxid, parent, level, name := parts[0], parts[1], parts[2], parts[4]
		// if !strings.HasPrefix(level, "S") {
		// 	continue
		// }
		if _, ok := m[taxid]; ok {
			return nil, fmt.Errorf("duplicate taxid: %q", taxid)
		}
		m[taxid] = &taxEntry{Name: name, Level: level, ParentTID: parent}
	}

	for tid, e := range m {
		if e.Level == "R" { // Root has no parents.
			continue
		}
		p := m[e.ParentTID]
		p.ChildTIDs = append(p.ChildTIDs, tid)
	}

	return m, nil
}
