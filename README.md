# SPM.jl

## Introduction

`SPM.jl` is a [Julia](https://julialang.org/) package alongside
the [SPM software](https://www.fil.ion.ucl.ac.uk/spm/). It has been
designed for the analysis of brain imaging data sequences.

## Installation

`SPM.jl` is available for Julia v1.0 and later releases, and
requires the packages `SpecialFunctions` and `NIfTI`.

```julia
import Pkg;
Pkg.add("SpecialFunctions")
Pkg.add("NIfTI")
```

It requires a shared library `libSPM` (see
SPM [Makefile](https://github.com/spm/spm/blob/main/src/Makefile)).
A compiled library `libSPM.so` is provided for 64-bit Linux environments.

## Usage

Example of between modality coregistration using information theory
(see [spm_coreg.m](https://github.com/spm/spm/blob/main/spm_coreg.m)).

If needed, you can download test data [here](https://github.com/spm/spm-notebooks/tree/main/data).

```julia
julia> include("Coreg.jl");

julia> P1 = "./T1w.nii";

julia> P2 = "./bold.nii";

julia> x,o = Coreg.run(P1,P2,[4],"nmi")
bold.nii -> T1w.nii w. nmi, 4mm samp. & 5x5 hist. smo.
0.4885   0        0        0        0        0        | -1.0011
0.4885   0.5516   0        0        0        0        | -1.0011
0.4885   0.5516   0        0        0        0        | -1.0011
0.4885   0.5516   0        0.1212   0        0        | -1.0012
0.4885   0.5516   0        0.1212   0        0        | -1.0012
0.4885   0.5516   0        0.1212   0        0        | -1.0012
0.497    0.5516   0        0.1212   0        0        | -1.0012
0.497    0.8482   0        0.1212   0        0        | -1.0012
0.497    0.8482   0.03426  0.1212   0        0        | -1.0012
0.497    0.8482   0.03426  0.1219   0        0        | -1.0012
0.497    0.8482   0.03426  0.1219   0        0        | -1.0012
0.497    0.8482   0.03426  0.1219   0        0        | -1.0012
0.5084   1.244    0.08003  0.1227   0        0        | -1.0012
0.5084   1.244    0.08003  0.1227   0        0        | -1.0012
0.5085   1.246    0.08019  0.1227   0        0        | -1.0012
0.5085   1.246    0.08019  0.1227   0        0        | -1.0012
0.5085   1.246    0.08019  0.1227   0        0        | -1.0012
0.5085   1.246    0.08019  0.1227   0        0        | -1.0012
0.5085   1.246    0.08019  0.1227   0        0        | -1.0012
([0.5084690641189673, 1.2458319387515466, 0.08018609030614379, 0.12270502915244379, 0.0, 0.0], -1.0012251f0)

```
