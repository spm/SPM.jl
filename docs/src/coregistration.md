
## Load package
```@example coreg
using SPM, NIfTI, CairoMakie, Tar, CodecZlib
```

## Download NIfTI and list path

!!! Note
    Artifact system is not working, we will download that by hand


```@example coreg
path_tar = download("https://github.com/spm/spm-notebooks/archive/2a09e1845ad6ec2ec6b9d0ca175cf8aa5d073e81.tar.gz",joinpath(@__DIR__,"spm_notebooks.tar.gz"))

tar_gz = open(path_tar)
tar = GzipDecompressorStream(tar_gz)
dir = Tar.extract(tar,joinpath(@__DIR__,"extracted"))
close(tar_gz)

notebooks_data = joinpath(dir,"spm-notebooks-2a09e1845ad6ec2ec6b9d0ca175cf8aa5d073e81","data")

const P1 = joinpath(notebooks_data, "T1w.nii")
const P2 = joinpath(notebooks_data, "bold.nii")
```

Now we can load them with the NIfTI package :
```@example coreg
T1 = niread(P1)
size(T1)
```

```@example coreg
bold = niread(P2)
size(bold)
```

```@example coreg
f=Figure()
ax=Axis(f[1,1],title="T1 slice 68", aspect = 1)
hidedecorations!(ax)
heatmap!(ax,abs.(T1.raw)[:,:,68])

ax=Axis(f[1,2],title="Bold slice 17", aspect = 1)
hidedecorations!(ax)
heatmap!(ax,abs.(bold.raw)[:,:,17])
f
```

## Run SPM

SPM works directly with nifti files path

```@example coreg
x,o = SPM.Coreg.run(P1,P2,[4],"nmi")
```


