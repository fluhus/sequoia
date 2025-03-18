package main

import (
	"testing"

	"golang.org/x/exp/slices"
)

func Test_ParseTaxIDs(t *testing.T) {
	tests := []struct {
		arg     string
		want    []int
		wantErr bool
	}{
		{"", nil, true},
		{"[111]", []int{111}, false},
		{"[22, 333]", []int{22, 333}, false},
		{"[4, 55, 6666]", []int{4, 55, 6666}, false},
		{"[4", nil, true},
		{"6666]", nil, true},
	}
	for _, test := range tests {
		got, err := Entry{}.ParseTaxIDs(test.arg)
		if (err != nil) != test.wantErr {
			t.Errorf("parseTaxIDs() error = %v, wantErr %v", err, test.wantErr)
			return
		}
		if !slices.Equal(got, test.want) {
			t.Errorf("parseTaxIDs(%q) = %v, want %v", test.arg, got, test.want)
		}
	}
}
