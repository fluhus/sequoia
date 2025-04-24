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
)

var (
	rgx    = flagx.Regexp("x", nil, "Pattern to search for")
	names  = flag.String("n", "", "Comma-separated names to search for (case-insensitive)")
	outDir = flag.String("o", "", "Output directory")
	inFile = flagx.FileExists("i", "", "Input fasta to search in")
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
		if !(*rgx).Match(fa.Name) {
			continue
		}
		n++
		txt, _ := fa.MarshalText()
		name := prefixRE.FindString(string(fa.Name))
		lens[name] = len(fa.Sequence)
		common.Die(os.WriteFile(*outDir+"/"+name+".fasta", txt, 0o644))
	}
	common.Die(jio.Write(*outDir+"/lens.json", lens))
	pt.Done()
}
