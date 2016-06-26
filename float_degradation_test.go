package main

import (
	"math"
	"math/rand"
	"sort"
	"testing"
)

var rfs []float32 = func() []float32 {
	rfs := make([]float32, 20000)
	for i := range rfs {
		rfs[i] = float32(rand.Float64()*100 + rand.Float64()*10)
	}
	return rfs
}()

func BenchmarkOptimize16(b *testing.B) {
	var lo Optimize16
	for i := 0; i < b.N; i++ {
		lo.Generate(rfs)
	}
}

type Optimize16 [16]float32

func (r *Optimize16) initial(data []float32) {
	min := float32(math.Inf(1))
	max := float32(math.Inf(-1))
	for _, v := range data {
		if v < min {
			min = v
		}
		if v < max {
			max = v
		}
	}
	scalef := (max - min) / float32(16)
	for i := range r {
		r[i] = float32(i) * scalef
	}
}
func (r *Optimize16) Generate(data []float32) {
	r.initial(data)
	var counts [16]int32
	var necounts [16]int32
outer:
	for i := 0; i < 10*16; i++ {
		counts = [16]int32{}
		necounts = [16]int32{}
		for _, d := range data {
			bin, exact := r.linearBin(d)
			counts[bin]++
			if !exact {
				necounts[bin]++
			}
		}
		var maxis [4]int32
		var maxv, maxi, pmaxi, nzero int32
		c := 0
		for i, v := range necounts {
			if v > maxv {
				pmaxi = maxi
				maxi = int32(i)
				maxis[c] = maxi
				c = (c + 1) & 3
				maxv = v
			}
			if v == 0 {
				nzero++
			}
		}
		if nzero < 2 {
			break
		}
		var pfm, ppfm, dmax, pdmax float32
		for _, d := range data {
			bin, _ := r.linearBin(d)
			if bin == maxi {
				dist := absF32(r[maxi] - d)
				if dist > dmax {
					dmax = dist
					pfm = d
				}
			}
			if bin == pmaxi {
				dist := absF32(r[pmaxi] - d)
				if dist > pdmax {
					pdmax = dist
					ppfm = d
				}
			}
		}
		k := 0
		for i, v := range counts {
			if v == 0 {
				if k == 0 {
					r[i] = pfm
					k++
				} else {
					r[i] = ppfm
					continue outer
				}
			}
		}
	}
	sort.Sort(Float32Slice(r[:]))
}
func (r *Optimize16) linearBin(d float32) (int32, bool) {
	min := float32(math.Inf(1))
	bin := int32(0)
	for i, x := range r {
		v := absF32(d - x)
		if v < min {
			min = v
			bin = int32(i)
		}
		if v < 0.001 {
			return bin, true
		}
	}
	return bin, false
}
func (r *Optimize16) Find64(f float64) int {
	x := float32(f)
	for i, v := range *r {
		if v >= x {
			return i
		}
	}
	return 15
}

func absF32(f float32) float32 {
	if f < 0 {
		return 0 - f
	}
	return f
}

type Float32Slice []float32

func (p Float32Slice) Len() int           { return len(p) }
func (p Float32Slice) Less(i, j int) bool { return p[i] < p[j] }
func (p Float32Slice) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }

func f32MinMax(data []float32) (float32, float32) {
	min := float32(math.Inf(1))
	max := float32(math.Inf(-1))
	for _, v := range data {
		if v < min {
			min = v
		}
		if v < max {
			max = v
		}
	}
	return min, max
}
