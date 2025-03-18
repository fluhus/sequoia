package main

import (
	"encoding/csv"
	"fmt"
	"lab/common"
	"os"

	"github.com/fluhus/gostuff/csvdec"
	"github.com/fluhus/gostuff/jio"
)

func main() {
	mapFile, inFile, outFile := os.Args[1], os.Args[2], os.Args[3]

	mm := map[string][]string{}
	common.Die(jio.Read(mapFile, &mm))
	m := map[string]string{}
	for k, v := range mm {
		for _, x := range v {
			m[x] = k
		}
	}
	fmt.Println(len(m))

	abnd := map[string]float64{}
	for l, err := range csvdec.File[abndLine](inFile, toTSV) {
		common.Die(err)
		abnd[l.Name] = l.Abnd
	}
	abnd2 := map[string]float64{}
	for k, v := range abnd {
		abnd2[m[k]] += v
	}
	fmt.Println(len(abnd), len(abnd2))

	common.Die(jio.Write(outFile, abnd2))
}

type abndLine struct {
	Name string
	Abnd float64
}

func toTSV(r *csv.Reader) {
	r.Comma = '\t'
}
