// Package common provides common utilities.
package common

import (
	"fmt"
	"math"
	"os"

	"github.com/fluhus/gostuff/gnum"
	"github.com/fluhus/gostuff/sets"
	"github.com/fluhus/gostuff/snm"
	"golang.org/x/exp/constraints"
	"golang.org/x/exp/maps"
)

// Die prints the error and exits if the error is non-nil.
func Die(err error) {
	if err != nil {
		fmt.Fprintln(os.Stderr, "ERROR:", err)
		os.Exit(2)
	}
}

// NTiles returns n+1 quantiles of s.
func NTiles[S ~[]E, E constraints.Ordered](n int, s S) S {
	result := make(S, n+1)
	for i := 0; i <= n; i++ {
		j := int(math.Round(float64(i) / float64(n) * float64(len(s)-1)))
		result[i] = s[j]
	}
	return result
}

func Perc[A, B gnum.Number](a A, b B, point int) string {
	if point < 0 {
		return fmt.Sprint(100*float64(a)/float64(b), "%")
	}
	return fmt.Sprintf("%.*f%%", point, 100*float64(a)/float64(b))
}

func NewDefaultMapPtr[K comparable, V any]() snm.DefaultMap[K, *V] {
	return snm.NewDefaultMap(func(k K) *V {
		var v V
		return &v
	})
}

func Dedup[T comparable](a []T) []T {
	if len(a) <= 1 {
		return a
	}
	return maps.Keys(sets.Of(a...))
}

func DedupSlice[T comparable](a [][]T) {
	for i := range a {
		a[i] = Dedup(a[i])
	}
}

func DedupValues[K comparable, V comparable](m map[K][]V) {
	for k, v := range m {
		m[k] = Dedup(v)
	}
}
