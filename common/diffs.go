package common

import "golang.org/x/exp/constraints"

func ToDiffs[I constraints.Integer](a []I) {
	last := len(a) - 1
	if last == -1 {
		return
	}
	for i := range a[1:] {
		a[last-i] -= a[last-i-1]
	}
}

func FromDiffs[I constraints.Integer](a []I) {
	if len(a) == 0 {
		return
	}
	for i := range a[1:] {
		a[i+1] += a[i]
	}
}
