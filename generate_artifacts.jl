using ArtifactUtils, Artifacts

add_artifact!(
    joinpath(@__DIR__, "Artifacts.toml"),
    "spm_notebooks",
    "https://github.com/spm/spm-notebooks/archive/2a09e1845ad6ec2ec6b9d0ca175cf8aa5d073e81.tar.gz",
    force=true,
    lazy=true,
)
