import SPM
using LazyArtifacts, Artifacts

const notebooks_data = joinpath(artifact"spm_notebooks", "spm-notebooks-2a09e1845ad6ec2ec6b9d0ca175cf8aa5d073e81", "data")

const P1 = joinpath(notebooks_data, "T1w.nii")
const P2 = joinpath(notebooks_data, "bold.nii")

x,o = SPM.Coreg.run(P1,P2,[4],"nmi")
A   = SPM.Coreg.spm_matrix(x)
