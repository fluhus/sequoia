// Package luna provides sample metadata handling functions.
package luna

import (
	"fmt"
	"path/filepath"
	"regexp"
	"strings"
)

// SampleNameMapping maps from v2 sample prefix to standardized sample name.
var SampleNameMapping = map[string]string{
	"A1": "Euro_Tur_111622",
	"A2": "Inh_Tur_111622",
	"B1": "Euro_Tur_113022",
	"B2": "Inh_Tur_113022",
	"C1": "Euro_Tur_040722",
	"C2": "Inh_Tur_040722",
	"D1": "Euro_Wod_111622",
	"D2": "Inh_Wod_111622",
	"E1": "Euro_Wod_113022",
	"H1": "Inh_Wod_113022",
	"F1": "Euro_Wod_041522",
	"G1": "Inh_Wod_041522",
	"H2": "Euro_LB_111622",
	"E2": "Inh_LB_111622",
	"A3": "Euro_LB_113022",
	"F2": "Inh_LB_113022",
	"B3": "Euro_LB_041422",
	"G2": "Inh_LB_041422",
}

// GroupNameMapping maps from group name to its name in the publication.
var GroupNameMapping = map[string]string{
	"1.Euro": "v1 solid",
	"1.Inh":  "v1 influent",
	"2.Euro": "v2 solid",
	"2.Inh":  "v2 influent",
}

// SampleOrdering is the list of samples, by run order.
var SampleOrdering = []string{
	"Inh_Tur_113022",
	"Inh_Wod_111622",
	"Inh_Tur_040722",
	"Inh_Wod_113022",
	"Inh_LB_113022",
	"Euro_Tur_111622",
	"Euro_Tur_113022",
	"Inh_Tur_111622",
	"Inh_LB_041422",
	"Inh_LB_111622",
	"Euro_LB_111622",
	"Euro_Wod_113022",
	"Euro_LB_113022",
	"Euro_Tur_040722",
	"Euro_LB_041422",
	"Euro_Wod_111622",
	"Euro_Wod_041522",
	"Inh_Wod_041522",
	"A1_S1_L008",
	"A2_S9_L008",
	"A3_S17_L008",
	"B1_S2_L008",
	"B2_S10_L008",
	"B3_S18_L008",
	"C1_S3_L008",
	"C2_S11_L008",
	"D1_S4_L008",
	"D2_S12_L008",
	"E1_S5_L008",
	"E2_S13_L008",
	"F1_S6_L008",
	"F2_S14_L008",
	"G1_S7_L008",
	"G2_S15_L008",
	"H1_S8_L008",
	"H2_S16_L008",
	"Undetermined_S0_L008",
}

// FixName converts a file name to standardized sample name.
func FixName(name string) string {
	s := filepath.Base(name)
	for _, p := range []string{".json", ".tax", ".krk", ".brk", ".brkraw", ".nreads"} {
		s = strings.TrimSuffix(s, p)
	}
	if strings.HasPrefix(s, "Euro_") || strings.HasPrefix(s, "Inh_") {
		return "1." + s
	}
	suf := SampleNameMapping[s[:2]]
	if suf == "" {
		panic(fmt.Sprintf("bad name: %q", name))
	}
	return "2." + SampleNameMapping[s[:2]]
}

// Detects the group part in a sample name.
var groupRE *regexp.Regexp

// Populates groupRE.
func init() {
	str := &strings.Builder{}
	str.WriteByte('^')
	for k := range GroupNameMapping {
		if str.Len() > 1 {
			str.WriteByte('|')
		}
		str.WriteString(regexp.QuoteMeta(k))
	}
	groupRE = regexp.MustCompile(str.String())
}

// SampleGroup extracts group name from a sample name.
func SampleGroup(s string) string {
	g := groupRE.FindString(s)
	if g == "" {
		panic(fmt.Sprintf("bad sample name: %q", s))
	}
	return g
}
