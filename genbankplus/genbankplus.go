// Package genbankplus makes some information in GenBank
// easily accessible.
package genbankplus

import (
	"fmt"
	"iter"
	"regexp"
	"strings"

	"github.com/fluhus/biostuff/formats/genbank"
)

type GenBankPlus struct {
	genbank.GenBank
	TaxID   string
	Host    string
	MolType string
}

func Iter(it iter.Seq2[*genbank.GenBank, error],
) iter.Seq2[*GenBankPlus, error] {
	return func(yield func(*GenBankPlus, error) bool) {
		for gb, err := range it {
			if err != nil {
				if !yield(&GenBankPlus{GenBank: *gb}, err) {
					return
				}
				continue
			}
			if !yield(ToPlus(gb)) {
				return
			}
		}
	}
}

func ToPlus(gb *genbank.GenBank) (*GenBankPlus, error) {
	p := &GenBankPlus{GenBank: *gb}

	// TaxID
	for _, f := range p.Features {
		taxID := f.Fields["db_xref"]
		if !strings.HasPrefix(taxID, "taxon:") {
			continue
		}
		taxID = taxID[6:] // Remove "taxon:".
		if taxID == "" {
			continue
		}
		if p.TaxID != "" && p.TaxID != taxID {
			p.TaxID = ""
			break
			// return p, fmt.Errorf("found two tax IDs: %q, %q", p.TaxID, taxID)
		}
		p.TaxID = taxID
	}

	// Host
	for _, f := range p.Features {
		host := f.Fields["host"]
		if host == "" {
			continue
		}
		if p.Host != "" && p.Host != host {
			p.Host = ""
			break
			// return p, fmt.Errorf("found two hosts: %q, %q", p.Host, host)
		}
		p.Host = host
	}

	// MolType
	molType := molTypeRE.FindStringSubmatch(gb.Locus)
	if molType == nil {
		return p, fmt.Errorf("could not find molecule-type in locus: %q",
			gb.Locus)
	}
	p.MolType = molType[1]

	return p, nil
}

// Captures the molecule-type part of the locus field.
var molTypeRE = regexp.MustCompile(`\s\d+ bp\s+(\S+)`)
