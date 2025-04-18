// Package genbankplus makes some information in GenBank
// easily accessible.
package genbankplus

import (
	"fmt"
	"iter"
	"regexp"
	"strconv"
	"strings"

	"github.com/fluhus/biostuff/formats/genbank"
)

type GenBankPlus struct {
	genbank.GenBank
	LocusName      string
	SequenceLength int
	MolType        string
	TaxID          string
	Host           string
}

func Iter(it iter.Seq2[*genbank.GenBank, error],
) iter.Seq2[*GenBankPlus, error] {
	return func(yield func(*GenBankPlus, error) bool) {
		for gb, err := range it {
			if err != nil {
				var gbp *GenBankPlus
				if gb != nil {
					gbp = &GenBankPlus{GenBank: *gb}
				}
				if !yield(gbp, err) {
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

	// LocusName, MolType
	locusParts := locusPartsRE.FindStringSubmatch(gb.Locus)
	if locusParts == nil {
		return p, fmt.Errorf("could not parse locus: %q", gb.Locus)
	}
	p.LocusName = locusParts[1]
	p.MolType = locusParts[3]

	// SequenceLength
	bp, err := strconv.Atoi(locusParts[2])
	if err != nil {
		return p, err
	}
	if bp < 1 {
		return p, fmt.Errorf("bad sequence length: %d", bp)
	}
	p.SequenceLength = bp

	return p, nil
}

// Captures the length and molecule-type parts of the locus field.
var locusPartsRE = regexp.MustCompile(`^(.+?)\s+(\d+) bp\s+(\S+)`)
