"""
    module SPMlib

This module accesses various compiled C functions used in the mex files of **SPM**.

"""
module SPMlib

export fmg3,cgs3,vel2mom,push,pull,identity,hist2,bound

const mwSize        = Csize_t
const mwSignedIndex = Csize_t
const libSPM        = joinpath(@__DIR__,"libSPM.so")

import LinearAlgebra:I

#import Libdl
#lib                = Libdl.dlopen(libSPM)
# This idea did not work because ccall is a keyword, rather than a function.
# The table is retained in case it us useful later on
const D = Dict(
         [(:get_bound,          (Cint,    ())),
          (:set_bound,          (Cvoid,   (Cint,))),
          (:fmg3_scratchsize,   (mwSize,  (Ptr{mwSize},Cint))),
          (:fmg3,               (Cvoid,   (Ptr{mwSize},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cdouble},Cint,Cint,Ptr{Cfloat},Ptr{Cfloat}))),
          (:cgs3,               (Cvoid,   (Ptr{mwSize},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cdouble},Cdouble,Cint,Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat}))),
          (:resize,             (Cvoid,   (Ptr{mwSize},Ptr{Cfloat},Ptr{mwSize},Ptr{Cfloat},Ptr{Cfloat}))),
          (:restrict_vol,       (Cvoid,   (Ptr{mwSize},Ptr{Cfloat},Ptr{mwSize},Ptr{Cfloat},Ptr{Cfloat}))),
          (:norm,               (Cdouble, (mwSize,Ptr{Cfloat}))),
          (:solve,              (Cvoid,   (Ptr{mwSize},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cfloat}))),
          (:LtLf,               (Cvoid,   (Ptr{mwSize},Ptr{Cfloat},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cfloat}))),
          (:fmg_scratchsize,    (Cint,    (Ptr{mwSize}))),
          (:trapprox,           (Cdouble, (Ptr{mwSize},Ptr{Cfloat},Ptr{Cdouble}))),
          (:vel2mom,            (Cvoid,   (Ptr{mwSize},Ptr{Cfloat},Ptr{Cdouble},Ptr{Cfloat}))),
          (:relax,              (Cvoid,   (Ptr{mwSize},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cdouble},Cint,Ptr{Cfloat}))),
          (:Atimesp,            (Cvoid,   (Ptr{mwSize},Ptr{Cfloat},Ptr{Cdouble},Ptr{Cfloat},Ptr{Cfloat}))),
          (:sumsq,              (Cdouble, (Ptr{mwSize},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cdouble},Ptr{Cfloat}))),
          (:kernel,             (Cvoid,   (Ptr{mwSize},Ptr{Cdouble},Ptr{Cfloat}))),
          (:composition,        (Cvoid,   (Ptr{mwSize},Ptr{mwSize},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat}))),
          (:composition_jacdet, (Cvoid,   (Ptr{mwSize},Ptr{mwSize},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat}))),
          (:push,               (Cvoid,   (Ptr{mwSize},mwSize,mwSize,Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat}))),
          (:pushc,              (Cvoid,   (Ptr{mwSize},mwSize,mwSize,Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat}))),
          (:pull,               (Cvoid,   (Ptr{mwSize},mwSize,mwSize,Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat}))),
          (:pullc,              (Cvoid,   (Ptr{mwSize},mwSize,mwSize,Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat}))),
          (:pushc_grads,        (Cvoid,   (Ptr{mwSize},Ptr{mwSize},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat}))),
          (:determinants,       (Cvoid,   (Ptr{mwSize},Ptr{Cfloat},Ptr{Cfloat}))),
          (:divergence,         (Cvoid,   (Ptr{mwSize},Ptr{Cfloat},Ptr{Cfloat}))),
          (:def2det,            (Cvoid,   (Ptr{mwSize},Ptr{Cfloat},Ptr{Cfloat},mwSignedIndex))),
          (:def2jac,            (Cvoid,   (Ptr{mwSize},Ptr{Cfloat},Ptr{Cfloat},mwSignedIndex))),
          (:invdef,             (Cvoid,   (Ptr{mwSize},Ptr{Cfloat},Ptr{mwSize},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat}))),
          (:grad,               (Cvoid,   (Ptr{mwSize},mwSize,Ptr{Cfloat},Ptr{Cfloat}))),
          (:hist2,              (Cvoid,   (Ptr{Cdouble},Ptr{Cuchar}, Ptr{Cuchar}, Ptr{mwSize}, Ptr{mwSize}, Ptr{Cdouble}, Ptr{Cfloat})))]);

## Removed because it did not work (ccall is a keyword, rather than a function)
#function ccall_args(s)
#    sym = Libdl.dlsym(lib, s)
#    return((sym,D[s][1],D[s][2]))
#end

"""
extern void pushc_grads(mwSize dmo[], mwSize dmy[], float def[], float J[], float pf[], float po[]);
extern void determinant(mwSize dm[], float J0[], float d[]);
extern void minmax_div(mwSize dm[], float v0[], double mnmx[]);
extern void divergence(mwSize dm[], float v0[], float dv[]);
extern void def2det(mwSize dm[], float *Y, float *J, mwSignedIndex s);
extern void def2jac(mwSize dm[], float *Y, float *J, mwSignedIndex s);
extern void invdef(mwSize dim_y[3], float y[], mwSize dim_iy[3], float iy[], float M1[4][3], float M2[4][3]);
extern void grad(mwSize dm[], mwSize N, float F[], float D[]);
extern void expm22(float a[], /*@out@*/ float l[]);
extern void expm33(float a[], /*@out@*/ float l[]);
extern void fmg3(mwSize n0[], float *a0, float *b0, double param[], int c, int nit,
                 float *u0, float *scratch);
extern void cgs3(mwSize dm[], float A[], float b[], double param[], double tol, int nit,
                 float x[], float r[], float p[], float Ap[]);
extern void resize(mwSize na[], float *a, mwSize nc[], float *c, float *b);
extern void restrict_vol(mwSize na[], float *a, mwSize nc[], float *c, float *b);
extern double norm(mwSize m, float a[]);
extern mwSize fmg3_scratchsize(mwSize n0[], int use_hessian);
extern void fmg(mwSize n0[], float *a0, float *b0, double param[], double scal[], int c, int nit, float *u0, float *scratch);
extern void solve(mwSize dm[], float a[], float b[], double s[], double scal[], float u[]);
extern void LtLf(mwSize dm[], float f[], double s[], double scal[], float g[]);
extern int fmg_scratchsize(mwSize n0[]);
extern double trapprox(mwSize dm[], float a[], double s[]);
extern void vel2mom(mwSize dm[], float f[], double s[], float g[]);
extern void relax(mwSize dm[], /*@null@*/ float a[], float b[], double s[], int nit, float u[]);
extern void Atimesp(mwSize dm[], /*@null@*/ float A[], double param[], float p[], float Ap[]);
extern double sumsq(mwSize dm[], float a[], float b[], double s[], float u[]);
extern void kernel(mwSize dm[], double s[], float f[]);
"""



function bound()
    return(ccall((:get_bound,libSPM),Cint,()))
end

function bound(b::Int)
    if b<0 || b>1
        error("inappropriate boundary condition")
    end
    ccall((:set_bound,libSPM),Cvoid,(Cint,),b)
    return(b)
end # bound

function fix_lambda(lambda::Array{Float64})
    if prod(size(lambda))!=3+5
        error("incorrect number of elements in lambda")
    end
    lambda = Float64.(lambda)
end # bound


"""
    fmg3(H,g,lambda,c,nit)

Solve equations using the Full Multigrid method.

# Example
```
H            = zeros(Float32,100,100,100,6);
H[:,:,:,1:3] = exp.(randn(Float32,100,100,100,3));
g            = randn(Float32,100,100,100,3);
lambda       = [1 1 1  0 0.1 1 0 0];
u            = SPMlib.fmg3(H,g,lambda,5,3);
```
"""
function fmg3(H::Array{Float32,4},g::Array{Float32,4},lambda::Array{Float64},c=2::Integer,nit=2::Integer)

    lambda = fix_lambda(lambda)

    # Check dimensions match
    da = size(H)
    db = size(g)
    if !(db[1:3]==da[1:3] && da[4]==6 && db[4]==3)
        error("dimension mismatch")
    end

    # Allocate scratch memory
    n0      = [db[1:3]...]

    s       = ccall((:fmg3_scratchsize,libSPM),Cintmax_t,(Ptr{Cintmax_t},Cint),
                                                          n0,            Int32(1))
    scratch = zeros(Float32,s)

    # Run fmg3 to get u as output
    u       = zeros(Float32,db)

    ccall((:fmg3,libSPM),Cvoid,
        (Ptr{Csize_t}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cdouble}, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}),
        n0,            H,           g,           lambda,       c,    nit,  u,           scratch)
    return(u)
end # fmg3

"""
    cgs3(H,g,lambda,tol,nit)

Solve equations using the conjugate gradient method.

# Example
```
H            = zeros(Float32,100,100,100,6);
H[:,:,:,1:3] = exp.(randn(Float32,100,100,100,3));
g            = randn(Float32,100,100,100,3);
lambda       = [1 1 1  0 0.1 1 0 0];
u            = SPMlib.cgs3(H,g,lambda,1e-4,1000);
```
"""
function cgs3(H::Array{Float32,4},g::Array{Float32,4},lambda::Array{Float64},tol=1e-4::AbstractFloat,nit=128::Integer)

    lambda = fix_lambda(lambda)

    # Check dimensions match
    da = size(H)
    db = size(g)
    if !(db[1:3]==da[1:3] && da[4]==6 && db[4]==3)
        error("dimension mismatch")
    end
    v0       = zeros(Float32,db)

    # Allocate scratch memory
    scratch1 = zeros(Float32,prod(db))
    scratch2 = zeros(Float32,prod(db))
    scratch3 = zeros(Float32,prod(db))

    n0      = [db[1:3]...]

    ccall((:cgs3,libSPM), Cvoid,
        (Ptr{mwSize},  Ptr{Cfloat},Ptr{Cfloat},Ptr{Cdouble},Cdouble,Cint,Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat},Ptr{Cfloat}),
         [db[1:3]...], H,          g,          lambda,      tol,    nit, v0,         scratch1,   scratch2,   scratch3)
    return(v0)
end # cgs3


function vel2mom(v0::Array{Float32,4},lambda::Array{Float64})
    lambda = fix_lambda(lambda)
    dm     = size(v0)
    u0     = zeros(Float32,dm)
    ccall((:vel2mom,libSPM), Cvoid, (Ptr{mwSize},Ptr{Cfloat},Ptr{Cdouble},Ptr{Cfloat}), [dm...],v0,lambda,u0);
    return(u0)
end # vel2mom


function pull(f0::Array{Float32},psi::Array{Float32,4},m=0)
    dpsi = size(psi);
    if dpsi[4]!=3
        error("wrong number of components")
    end
    df   = size(f0);
    dm0  = [df[1:3]...];
    m1   = prod(dpsi[1:3]);
    n    = prod(df[4:end]); # Should not work in most cases - but it does
    dm1  = (dpsi[1:3]..., df[4:end]...);
    f1   = zeros(Float32,dm1);
    code = (m!=0 ? 1 : 3)

    ccall((:pushpull,libSPM),Cvoid,
        (Ptr{Csize_t}, Csize_t, Csize_t, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Nothing}, Ptr{Cfloat}, Cuint),
        dm0,           m1,      n,       psi,         f0,          C_NULL,       f1,          code);
    return(f1)
end # pull


function push(f1::Array{Float32},psi::Array{Float32,4},dm1::NTuple{3,Integer},m=0)
    dpsi = size(psi);
    if dpsi[4]!=3 error("wrong number of components") end
    df   = size(f1);
    if ~all(df[1:3]==dpsi[1:3]) error("unmatched dimensions") end
    dm0  = [df[1:3]...];
    m1   = prod(dpsi[1:3]);
    n    = prod(df[4:end]); # Should not work in most cases - but it does
    f0   = zeros(Float32,dm1);
    code = (m!=0 ? 0 : 2)

    ccall((:pushpull,libSPM), Cvoid,
        (Ptr{Csize_t}, Csize_t, Csize_t, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Nothing}, Ptr{Cfloat}, Cuint),
        dm0,           m1,      n,       psi,         f0,          C_NULL,       f1,          code);
    return(f0)
end # push


"""
    H = hist2(f0, f1, M, s)

2D joint histogram for mutual information affine registration.

Inputs:
  * ```f0::Array{UInt8,3}```                        - Fixed 3D image
  * ```f1::Array{UInt8,3}```                        - Moved 3D image
  * ```M::Matrix{Float64}=Matrix{Float64}(I,4,4)``` - Affine mapping between images (4 × 4)
  * ```sep::Array{Float64}=[1.0, 1.0, 1.0]```       - Approximate separation between samples (voxels)

Output:
  * ```H::Array{Float64,2}```                       - Joint intensity histogram (256 × 256)

"""
function hist2(f0::Array{UInt8,3}, f1::Array{UInt8,3}, M::Matrix{Float64}=Matrix{Float64}(I,4,4), sep::Array{Float64}=[1.0, 1.0, 1.0])
    dm0 = size(f0)
    dm1 = size(f1)
    if length(sep) ==1
        sep = [1.0, 1.0, 1.0]*sep
    else
        sep = [sep[1:min(end,3)]' ones(Float32,1,max(3-length(sep),0))]
    end
    H = zeros(Float64,(256,256))
    ccall((:hist2,libSPM), Cvoid, (Ptr{Cdouble},Ptr{Cuchar}, Ptr{Cuchar}, Ptr{mwSize}, Ptr{mwSize}, Ptr{Cdouble}, Ptr{Cfloat}),
                                   Float64.(M), f0,          f1,          [dm0...],    [dm1...],    H,            Float32.(sep));
    return(H)
end # hist2


function grad(f0::Array{Float32})
    dm  = size(f0)
    dm  = (length(dm)<3 ? (dm...,ones(Integer,3-length(dm))...) : dm)
    dmi = [dm[1:3]...]
    n   = prod(dm[4:end])
    dmo = (dm[1:3]...,3,dm[4:end]...)
    gf  = zeros(Float32,dmo)
    ccall((:grad,libSPM),Cvoid,(Ptr{mwSize},mwSize,Ptr{Cfloat},Ptr{Cfloat}), dmi,n,f0,gf)
    return(gf)
end # grad


function identity(d::NTuple{3,Int64})
    psi = cat(reshape(Float32.(collect(1:d[1])),(d[1],1,1)).+zeros(Float32,(1,d[2],d[3])),
              reshape(Float32.(collect(1:d[2])),(1,d[2],1)).+zeros(Float32,(d[1],1,d[3])),
              reshape(Float32.(collect(1:d[3])),(1,1,d[3])).+zeros(Float32,(d[1],d[2],1)),dims=4);
    return(psi)
end # identity

end # module

