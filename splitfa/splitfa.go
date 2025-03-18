// Splits fasta files.
package main

import (
	"flag"
	"fmt"
	"lab/common"
	"lab/lazy"
	"runtime/debug"
	"strings"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/flagx"
	"github.com/fluhus/gostuff/ptimer"
	"github.com/fluhus/gostuff/snm"
)

var (
	inFile  = flag.String("i", "", "Input fasta file")
	outFile = flag.String("o", "", "Output file with % for the index or group")
	nseqs   = flag.Int("n", 0, "Number of sequences per file")
	groupRE = flagx.Regexp("g", nil, "Regular expression by which to split")
	ignore  = flag.Bool("ignore", false, "Ignore names that don't match regex")
	appnd   = flag.Bool("append", false, "Append to existing files instead of overwrite")
)

func main() {
	flag.Parse()

	if *nseqs != 0 && *groupRE != nil {
		common.Die(fmt.Errorf("only one of -n and -g is allowed"))
	}
	if *nseqs == 0 && *groupRE == nil {
		common.Die(fmt.Errorf("either -n or -g must be set"))
	}

	debug.SetGCPercent(20)

	if *nseqs != 0 {
		var fas []*fasta.Fasta
		i := 1
		pt := ptimer.New()
		for fa, err := range fasta.File(*inFile) {
			pt.Inc()
			common.Die(err)
			fas = append(fas, fa)
			if len(fas) >= *nseqs {
				fout, err := aio.Create(strings.Replace(*outFile, "%", fmt.Sprint(i), 1))
				common.Die(err)
				for _, fa := range fas {
					common.Die(fa.Write(fout))
				}
				fout.Close()
				i++
				fas = fas[:0]
			}
		}
		if len(fas) > 0 {
			fout, err := aio.Create(strings.Replace(*outFile, "%", fmt.Sprint(i), 1))
			common.Die(err)
			for _, fa := range fas {
				common.Die(fa.Write(fout))
			}
			fout.Close()
		}
		pt.Done()
	} else { // Group regex.
		re := *groupRE
		fout := snm.NewDefaultMap(func(s string) *lazy.Writer {
			if *appnd {
				return lazy.Append(s, 1<<17 /* 128KB */)
			} else {
				return lazy.Create(s, 1<<17 /* 128KB */)
			}
		})
		pt := ptimer.NewFunc(func(i int) string {
			return fmt.Sprintf("%v (%v files)", i, len(fout.M))
		})
		for fa, err := range fasta.File(*inFile) {
			pt.Inc()
			common.Die(err)
			g := re.FindStringSubmatch(string(fa.Name))
			if len(g) == 0 {
				if *ignore {
					continue
				}
				common.Die(fmt.Errorf("name did not match regex: %q", fa.Name))
			}
			m := g[0]
			if len(g) > 1 {
				m = g[1]
			}
			f := strings.Replace(*outFile, "%", m, 1)
			common.Die(fa.Write(fout.Get(f)))
		}
		for _, v := range fout.M {
			common.Die(v.Flush())
		}
		pt.Done()
	}
}
