// Extracts viral sequences from refseq and outputs them
// as single files, for coverage analysis.
package main

import (
	"flag"
	"fmt"
	"lab/common"
	"os"
	"regexp"
	"strings"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/gostuff/flagx"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/ptimer"
	"github.com/fluhus/gostuff/sets"
)

var (
	rgx    = flagx.Regexp("x", nil, "Pattern to search for")
	names  = flag.String("n", "", "Comma-separated names to search for (case-insensitive)")
	outDir = flag.String("o", "", "Output directory")
	inFile = flagx.FileExists("i", "", "Input fasta to search in")
	jsFile = flagx.FileExists("j", "", "Input JSON with names")
)

func main() {
	flag.Parse()
	if *names != "" {
		parts := strings.Split(*names, ",")
		fmt.Printf("Searching for: %q\n", parts)
		for i := range parts {
			parts[i] = regexp.QuoteMeta(parts[i])
		}
		pat := "(?i:" + strings.Join(parts, "|") + ")"
		*rgx = regexp.MustCompile(pat)
	}
	wl, err := readListFromJSON()
	common.Die(err)

	common.Die(os.MkdirAll(*outDir, 0o755))

	n := 0
	pt := ptimer.NewFunc(func(i int) string {
		return fmt.Sprintf("%d (%d)", i, n)
	})
	lens := map[string]int{}
	prefixRE := regexp.MustCompile(`^\S+`)
	for fa, err := range fasta.File(*inFile) {
		common.Die(err)
		pt.Inc()
		name := prefixRE.FindString(string(fa.Name))
		if !wl.Has(name) && !(*rgx).Match(fa.Name) {
			continue
		}
		n++
		txt, _ := fa.MarshalText()
		lens[name] = len(fa.Sequence)
		common.Die(os.WriteFile(*outDir+"/"+name+".fasta", txt, 0o644))
	}
	common.Die(jio.Write(*outDir+"/lens.json", lens))
	pt.Done()

	if n == 0 {
		common.Die(fmt.Errorf("no matching species found"))
	}
}

func readListFromJSON() (sets.Set[string], error) {
	s := sets.Set[string]{}
	if *jsFile == "" {
		return s, nil
	}
	var items []jItem
	if err := jio.Read(*jsFile, &items); err != nil {
		return nil, err
	}
	for _, x := range items {
		if (*rgx).MatchString(x.Source) || (*rgx).MatchString(x.Organism) {
			s.Add(x.Version)
		}
	}
	return s, nil
}

type jItem struct {
	Source   string
	Organism string
	Version  string
}
