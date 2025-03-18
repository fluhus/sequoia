package main

import (
	"math/rand/v2"
	"testing"

	"github.com/fluhus/gostuff/snm"
)

func BenchmarkShuffle(b *testing.B) {
	nums := snm.Slice(10000, func(i int) int { return i })
	b.Run("my", func(b *testing.B) {
		for range b.N {
			snm.Shuffle(nums)
		}
	})
	b.Run("std", func(b *testing.B) {
		for range b.N {
			rand.Shuffle(len(nums), func(i, j int) {
				nums[i], nums[j] = nums[j], nums[i]
			})
		}
	})
}
