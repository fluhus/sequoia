// Parses NCBI's XML dump.
package main

import (
	"encoding/xml"
	"fmt"
	"lab/common"
	"lab/config"
	"runtime/debug"
	"strings"
	"time"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/jio"
)

const (
	dir = config.DataDir + "/ncbi/"
)

func main() {
	debug.SetGCPercent(20)

	fmt.Println("Parsing")
	t := time.Now()
	var x struct {
		Items []struct {
			ID           string `xml:"INSDSeq_primary-accession"`
			Name         string `xml:"INSDSeq_organism"`
			Tax          string `xml:"INSDSeq_taxonomy"`
			Seq          string `xml:"INSDSeq_sequence"`
			MolType      string `xml:"INSDSeq_moltype"`
			FeatureTable struct {
				Features []struct {
					QualifierList struct {
						Qualifiers []struct {
							Name  string `xml:"INSDQualifier_name"`
							Value string `xml:"INSDQualifier_value"`
						} `xml:"INSDQualifier"`
					} `xml:"INSDFeature_quals"`
				} `xml:"INSDFeature"`
			} `xml:"INSDSeq_feature-table"`
			Host string
		} `xml:"INSDSeq"`
	}

	f, err := aio.Open(dir + "sequence.gbc.xml.zst")
	common.Die(err)
	common.Die(xml.NewDecoder(f).Decode(&x))
	f.Close()
	fmt.Println(time.Since(t))

	fmt.Println(len(x.Items), "sequences")

	fmt.Println("Checking")
	t = time.Now()
	badSeq := 0
	withHost := 0
	for _, it := range x.Items {
		if !isNuc(it.Seq) {
			badSeq++
		}
	}
	for i, it := range x.Items {
		for _, f := range it.FeatureTable.Features {
			for _, q := range f.QualifierList.Qualifiers {
				if strings.ToLower(q.Name) == "host" {
					x.Items[i].Host = q.Value
				}
			}
		}
		x.Items[i].FeatureTable.Features = nil
		if x.Items[i].Host != "" {
			withHost++
		}
	}
	fmt.Println(time.Since(t))
	fmt.Println("Bad sequence:", badSeq)
	fmt.Println("With host:", withHost)

	if false {
		fmt.Println("Writing fasta")
		t = time.Now()
		f, err := aio.Create(dir + "sequence.fasta")
		common.Die(err)
		for _, it := range x.Items {
			if !isNuc(it.Seq) {
				continue
			}
			fa := &fasta.Fasta{Name: []byte(it.ID), Sequence: []byte(it.Seq)}
			common.Die(fa.Write(f))
		}
		f.Close()
		fmt.Println(time.Since(t))
	}

	fmt.Println("Writing metadata")
	t = time.Now()
	tax := map[string]map[string]string{}
	for _, it := range x.Items {
		m := map[string]string{
			"name": it.Name, "tax": it.Tax, "molType": it.MolType,
		}
		if it.Host != "" {
			m["host"] = it.Host
		}
		tax[it.ID] = m
	}
	if len(tax) != len(x.Items) {
		common.Die(fmt.Errorf("mismatching lengths: %d, %d",
			len(x.Items), len(tax)))
	}
	if true {
		common.Die(jio.Write(dir+"sequence.tax.json", tax))
	}
	fmt.Println(time.Since(t))
}

func isNuc(s string) bool {
	for _, x := range s {
		if !nucTable[x] {
			return false
		}
	}
	return true
}

var nucTable = make([]bool, 256)

func init() {
	nucTable['a'] = true
	nucTable['A'] = true
	nucTable['c'] = true
	nucTable['C'] = true
	nucTable['g'] = true
	nucTable['G'] = true
	nucTable['t'] = true
	nucTable['T'] = true
	nucTable['n'] = true
	nucTable['N'] = true
}
