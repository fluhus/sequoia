// Outputs lengths of fasta sequences.
package main

import (
	"fmt"
	"lab/common"
	"os"
	"path/filepath"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/gostuff/gnum"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/snm"
	"golang.org/x/exp/maps"
)

const (
	trimNames = true
)

func main() {
	if trimNames {
		fmt.Println("TRIMMING NAMES")
	}

	lens := map[string]int{}
	files, err := filepath.Glob(os.Args[1])
	common.Die(err)
	for _, file := range files {
		for fa, err := range fasta.File(file) {
			common.Die(err)
			if trimNames {
				for i, c := range fa.Name {
					if c == ' ' {
						fa.Name = fa.Name[:i]
						break
					}
				}
			}
			if lens[string(fa.Name)] != 0 {
				common.Die(fmt.Errorf("name appeared twice: %q", fa.Name))
			}
			lens[string(fa.Name)] = len(fa.Sequence)
		}
	}
	jio.Write(os.Args[2], lens)

	sorted := snm.Sorted(maps.Values(lens))
	fmt.Print(len(sorted), " [")
	const n = 20
	for i := range n + 1 {
		ii := gnum.Idiv((len(sorted)-1)*i, n)
		fmt.Print(sorted[ii], " ")
	}
	fmt.Println("]")
}
