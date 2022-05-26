#!/usr/bin/env julia
push!(LOAD_PATH, joinpath(pwd(),"..","src"))
import SPM

using Profile:@profile

P1  = download("https://github.com/spm/spm-notebooks/raw/main/data/T1w.nii")
P2  = download("https://github.com/spm/spm-notebooks/raw/main/data/bold.nii")

x,o = SPM.Coreg.run(P1,P2,[4],"nmi")
A   = SPM.Coreg.spm_matrix(x)
