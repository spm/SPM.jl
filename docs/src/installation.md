## Installation

`SPM.jl` is available for Julia v1.6 and later releases.  It can be installed
with [Julia built-in package
manager](https://julialang.github.io/Pkg.jl/stable/).  In a Julia session, after
entering the package manager mode with `]`, run the command

```julia
pkg> update
pkg> add https://github.com/spm/SPM.jl
```

To exit from the Pkg REPL mode press the backspace key, or `Ctrl + C`.

`SPM.jl` also requires a shared library `libSPM` (see
[SPM Makefile](https://github.com/spm/spm/blob/main/src/Makefile)).
A compiled library `libSPM.so` is provided for 64-bit Linux environments.