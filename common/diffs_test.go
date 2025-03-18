package common

import (
	"slices"
	"testing"
)

func TestToDiffs(t *testing.T) {
	tests := []struct {
		input []int
		want  []int
	}{
		{nil, nil},
		{[]int{3}, []int{3}},
		{[]int{1, 4, 9, 16, 25}, []int{1, 3, 5, 7, 9}},
	}
	for _, test := range tests {
		got := slices.Clone(test.input)
		ToDiffs(got)
		if !slices.Equal(got, test.want) {
			t.Errorf("ToDiffs(%v)=%v, want %v", test.input, got, test.want)
		}
	}
}

func TestFromDiffs(t *testing.T) {
	tests := []struct {
		input []int
		want  []int
	}{
		{nil, nil},
		{[]int{3}, []int{3}},
		{[]int{1, 3, 5, 7, 9}, []int{1, 4, 9, 16, 25}},
	}
	for _, test := range tests {
		got := slices.Clone(test.input)
		FromDiffs(got)
		if !slices.Equal(got, test.want) {
			t.Errorf("FromDiffs(%v)=%v, want %v", test.input, got, test.want)
		}
	}
}
