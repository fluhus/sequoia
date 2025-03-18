// Extracts individual reference genomes of enriched taxa.
package main

import (
	"fmt"
	"lab/common"

	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/sets"
)

func main() {
	tid2nc := map[string][]string{}
	common.Die(jio.Read("../../../Data/refseq/tid2nc.json", &tid2nc))
	// common.Die(jio.Read("species_krk.json", &tid2nc)) // Kraken's mapping doesn't have all of its output species.
	fmt.Println(len(tid2nc))

	enrichedAccs := sets.Set[string]{}
	{
		enrichedRaw := map[string][]string{}
		common.Die(jio.Read("species_vsp2.json", &enrichedRaw))
		for _, accs := range enrichedRaw {
			enrichedAccs.Add(accs...)
		}
	}
	fmt.Println(len(enrichedAccs))

	covSpecies := struct {
		Top  sets.Set[string]
		Vsp2 sets.Set[string]
	}{}
	common.Die(jio.Read("cov_species.json", &covSpecies))
	fmt.Println(len(covSpecies.Top), len(covSpecies.Vsp2),
		len(sets.Set[string]{}.AddSet(covSpecies.Top).AddSet(covSpecies.Vsp2)))

	mapping := map[string]string{} // NC to taxid
	foundTIDs := sets.Set[string]{}
	for tid, accs := range tid2nc {
		if !covSpecies.Top.Has(tid) && !covSpecies.Vsp2.Has(tid) {
			continue
		}
		for _, acc := range accs {
			if !covSpecies.Top.Has(tid) && !enrichedAccs.Has(acc) {
				continue
			}
			if mapping[acc] != "" {
				common.Die(fmt.Errorf("duplicate tid for acc %s: %s, %s",
					acc, mapping[acc], tid))
			}
			mapping[acc] = tid
			foundTIDs.Add(tid)
		}
	}
	fmt.Println(len(mapping))
	common.Die(jio.Write("taxref.json", mapping))

	for tid := range covSpecies.Top {
		if !foundTIDs.Has(tid) {
			fmt.Printf("Not found (top): %s\n", tid)
		}
	}
	for tid := range covSpecies.Vsp2 {
		if !foundTIDs.Has(tid) {
			fmt.Printf("Not found (vsp): %s\n", tid)
		}
	}
}
