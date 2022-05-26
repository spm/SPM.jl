#!julia 
include("Coreg.jl")
using Profile:@profile
P1="./test/data/T1w.nii"
P2="./test/data/bold.nii"
x,o = Coreg.run(P1,P2,[4],"nmi")

