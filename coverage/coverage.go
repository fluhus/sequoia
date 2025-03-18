// Creates coverage histogram from a SAM file.
package main

import (
	"fmt"
	"lab/common"
	"os"
	"path/filepath"

	"github.com/fluhus/biostuff/formats/sam"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/ptimer"
	"github.com/fluhus/gostuff/snm"
)

func main() {
	var lens map[string]int
	if len(os.Args) == 4 {
		lens = map[string]int{}
		// common.Die(jio.Read("../../../Data/refseq/viral.1.genomic.lens.json", &lens))
		common.Die(jio.Read(os.Args[3], &lens))
	}

	refToCov := snm.NewDefaultMap(func(s string) []int {
		if lens != nil {
			return make([]int, lens[s])
		} else {
			return nil
		}
	})
	type bad struct {
		want int
		got  int
	}
	var bads []bad
	files, _ := filepath.Glob(os.Args[1])
	fi := 0
	pt := ptimer.NewFunc(func(i int) string {
		return fmt.Sprintf("%v (file %v/%v)", i, fi, len(files))
	})
	for _, file := range files {
		fi++
		for sm, err := range sam.File(file) {
			common.Die(err)
			if lens != nil && lens[sm.Rname] == 0 {
				common.Die(fmt.Errorf("unrecognized reference: %q", sm.Rname))
			}
			pos, ln := sm.Pos-1, len(sm.Seq)
			a := refToCov.Get(sm.Rname)
			if lens != nil && lens[sm.Rname] <= pos {
				bads = append(bads, bad{want: lens[sm.Rname], got: pos})
			}
			for len(a) < pos+ln {
				a = append(a, 0)
			}
			for i := range ln {
				a[pos+i]++
			}
			refToCov.Set(sm.Rname, a)
			pt.Inc()
		}
	}
	pt.Done()
	fmt.Println(len(refToCov.M), "reference sequences")

	for _, b := range bads {
		fmt.Printf("W:%10d G:%10d\n", b.want, b.got)
	}

	jio.Write(os.Args[2], refToCov.M)
}
