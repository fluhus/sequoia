// Looks up uniref entries in the GO database.
package main

import (
	"encoding/csv"
	"fmt"
	"iter"
	"lab/common"
	"os"
	"regexp"
	"runtime/debug"
	"slices"
	"strings"

	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/gnum"
	"github.com/fluhus/gostuff/iterx"
	"github.com/fluhus/gostuff/ptimer"
	"golang.org/x/exp/maps"
)

const (
	goFile = "/dfs7/whitesonlab/alavon/Data/go/filtered_goa_uniprot_all.gaf.gz"
)

func main() {
	debug.SetGCPercent(20)

	inFile, outFile := os.Args[1], os.Args[2]

	m, err := loadQueries(inFile)
	common.Die(err)
	nEntries := gnum.Max(slices.Concat(maps.Values(m)...)) + 1
	fmt.Println(nEntries, "entries,", len(m), "IDs to search")

	found := make([][]string, nEntries)

	pt := ptimer.NewFunc(func(i int) string {
		n := 0
		for _, x := range found {
			if x != nil {
				n++
			}
		}
		return fmt.Sprintf("%v (found %v)", i, n)
	})

	for e, err := range iterEntries(goFile) {
		common.Die(err)
		if ii, ok := m[e.Acc]; ok {
			for _, i := range ii {
				found[i] = append(found[i], e.GoTerm)
			}
		}
		pt.Inc()
	}
	pt.Done()

	fout, err := aio.Create(outFile)
	common.Die(err)
	for _, x := range found {
		_, err := fmt.Fprintln(fout, strings.Join(x, ","))
		common.Die(err)
	}
	common.Die(fout.Close())
}

func loadQueries(file string) (map[string][]int, error) {
	re := regexp.MustCompile(`UniRef\d+_(\w+)`)
	m := map[string][]int{}
	i := -1
	for line, err := range iterx.LinesFile(file) {
		if err != nil {
			return nil, err
		}
		i++
		match := re.FindAllStringSubmatch(line, -1)
		for _, x := range match {
			m[x[1]] = append(m[x[1]], i)
		}
	}
	return m, nil
}

func iterEntries(file string) iter.Seq2[goEntry, error] {
	return func(yield func(goEntry, error) bool) {
		it := iterx.CSVFile(file, func(r *csv.Reader) {
			r.Comma = '\t'
			r.Comment = '!'
			r.FieldsPerRecord = -1
		})
		i := -1
		for line, err := range it {
			i++
			if err != nil {
				yield(goEntry{}, fmt.Errorf("line %v: %w", i, err))
				return
			}
			if len(line) < 5 {
				continue
			}
			e := goEntry{
				DB:       line[0],
				Acc:      line[1],
				Gene:     line[2],
				Relation: line[3],
				GoTerm:   line[4],
			}
			if !strings.HasPrefix(e.GoTerm, "GO:") {
				yield(goEntry{}, fmt.Errorf("line %v: bad GO term: %q",
					i, e.GoTerm))
				return
			}
			if !yield(e, nil) {
				return
			}
		}
	}
}

type goEntry struct {
	DB       string
	Acc      string
	Gene     string
	Relation string
	GoTerm   string
}
