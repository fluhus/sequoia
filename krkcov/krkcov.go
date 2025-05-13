// Filters output reads from kraken, keeping only ones that mapped
// to specific tax IDs.
package main

import (
	"fmt"
	"lab/common"
	"lab/config"
	"os"
	"regexp"
	"runtime/debug"

	"github.com/fluhus/biostuff/formats/fastq"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/sets"
	"github.com/fluhus/gostuff/snm"
)

const (
	krkFile = config.WSDataDir + "/krk_std2.json"
)

func main() {
	debug.SetGCPercent(20)

	taxids, err := childrenOf("2843396") // Jouyvirus.
	common.Die(err)
	fmt.Fprintln(os.Stderr, len(taxids), "taxids")

	re := regexp.MustCompile(`kraken:taxid\|(\d+)`)
	for fq, err := range fastq.Reader(os.Stdin) {
		common.Die(err)
		m := re.FindSubmatch(fq.Name)
		if m == nil {
			common.Die(fmt.Errorf("no kraken taxid: %q", fq.Name))
		}
		if taxids.Has(string(m[1])) {
			fq.Write(os.Stdout)
		}
	}
}

// Returns the children of the given tax ID,
// from the Kraken DB hierarchy.
func childrenOf(tid string) (sets.Set[string], error) {
	m := map[string]*taxEntry{}
	if err := jio.Read(krkFile, &m); err != nil {
		return nil, err
	}
	s := sets.Set[string]{}
	q := &snm.Queue[string]{}
	q.Enqueue(tid)
	for tid := range q.Seq() {
		x := m[tid]
		if x == nil {
			return nil, fmt.Errorf("tid not found: %q", tid)
		}
		s.Add(tid)
		// fmt.Println("Name:", x.Name)
		for _, ctid := range x.ChildTIDs {
			q.Enqueue(ctid)
		}
	}
	return s, nil
}

type taxEntry struct {
	Name      string
	Level     string
	ParentTID string
	ChildTIDs []string `json:",omitempty"`
	Accs      []string `json:",omitempty"`
}
