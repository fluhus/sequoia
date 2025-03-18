// Checks that an input SAM file has interleaved paired reads.
package main

import (
	"fmt"
	"lab/common"
	"os"

	"github.com/fluhus/biostuff/formats/fastq"
	"github.com/fluhus/gostuff/ptimer"
)

func main() {
	bad := 0
	pt := ptimer.NewFunc(func(i int) string { return fmt.Sprintf("%v (%v)", i, bad) })
	last := ""
	for fq, err := range fastq.File(os.Args[1]) {
		common.Die(err)
		if last == "" {
			last = string(fq.Name)
		} else {
			if last != string(fq.Name) {
				bad++
				last = string(fq.Name)
			} else {
				last = ""
			}
		}
		pt.Inc()
	}
	pt.Done()
	if bad > 0 {
		os.Exit(1)
	}
}
