#!julia 
include("Coreg.jl")
using Profile:@profile
P1="./IXI326-Guys-0907-T1.nii"
P2="./IXI326-Guys-0907-T2.nii"
x,o = Coreg.run(P1,P2,[4],"nmi")

