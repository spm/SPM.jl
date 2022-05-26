module Coreg
using SpecialFunctions:erf
using LinearAlgebra
using Random:seed!,rand,randn
using Printf

function convsame(x::Array{<:Real},c::Array{<:Real})
    dx  = size(x)
    dc  = size(c)
#   if length(dx)<length(dc)
#       dx = (dx...,ones(Int64,length(dc)-length(dx))...)
#   elseif length(dx)>length(dc)
#       dc = (dc...,ones(Int64,length(dx)-length(dc))...)
#   end
    o   = CartesianIndex(Int64.(floor.((dc.+1)./2))...)
    ret = zeros(Float32,dx)
    Ic  = CartesianIndices(dc)
    Ix  = CartesianIndices(dx)
    for ix in Ix
        for ic in Ic
            ir = ix+ic-o
            if in(ir,Ix) ret[ir] += x[ix] * c[ic] end
        end
    end
    return ret
end

function convn(x::Array{<:Real},c::Array{<:Real})
    dx  = size(x)
    dc  = size(c)
#   if length(dx)<length(dc)
#       dx = (dx...,ones(Int64,length(dc)-length(dx))...)
#   elseif length(dx)>length(dc)
#       dc = (dc...,ones(Int64,length(dx)-length(dc))...)
#   end
    o   = CartesianIndex((ones(Int64,length(dx))...))
    ret = zeros(Float32,dx.+dc.-1)
    Ic  = CartesianIndices(dc)
    Ix  = CartesianIndices(dx)
    for ix in Ix
        for ic in Ic
            ret[ix+ic-o] += x[ix] * c[ic]
        end
    end
    return ret
end

"""

    krn = smoothkern(fwhm,x,t)

Generate a Gaussian smoothing kernel

Inputs:
  * ```fwhm``` - full width at half maximum
  * ```x```    - position
  * ```t```    - either 0 (nearest neighbour) or 1 (linear).
           [Default: 1]
Returns:
  * ```krn```  - value of kernel at position x

---

For smoothing images, one should really convolve a Gaussian with a sinc
function. For smoothing histograms, the kernel should be a Gaussian
convolved with the histogram basis function used. This function returns
a Gaussian convolved with a triangular (1st degree B-spline) basis
function (by default). A Gaussian convolved with a hat function (0th
degree B-spline) can also be returned.

Note that the convolution kernel returned by this function differ from
the ones that other packages currently use for Gaussian smoothing -
particularly when the FWHM is small compared with the voxel dimensions.
The fact that SPM does it differently from other software does not mean
that it is wrong.

---
Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging
"""
function smoothkern(fwhm::Real,x,t::Integer=1)

    # Variance from FWHM
    s = (fwhm/sqrt(8*log(2)))^2+eps(Float64)

    # The simple way to do it. Not good for small FWHM
    # krn = (1/sqrt(2*pi*s)).*exp.(-(x.^2)./(2*s))

    if t==0
        # Gaussian convolved with 0th degree B-spline
        # int(exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t= -0.5..0.5)
        w1  = 1/sqrt(2.0*s)
        krn = 0.5.*(erf.(w1.*(x.+0.5)).-erf.(w1.*(x.-0.5)))
        krn[krn.<0.0] .= 0.0

    elseif t==1
        # Gaussian convolved with 1st degree B-spline
        #  int((1-t)*exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t= 0..1)
        # +int((t+1)*exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t=-1..0)
        w1  =  1/sqrt(2.0*s)
        w2  = -0.5/s
        w3  = sqrt(s/2.0/pi)
        krn = 0.5.*(erf.(w1.*(x.+1)).*(x.+1) .+ erf.(w1.*(x.-1)).*(x.-1) .- 2.0.*erf.(w1.*x   ).* x) .+
               w3.*(exp.(w2.*(x.+1).^2)      .+ exp.(w2.*(x.-1).^2)      .- 2.0.*exp.(w2.*x.^2))
        krn[krn.<0.0] .= 0.0

    else
        error("Only defined for nearest neighbour and linear interpolation.");
        # This should probably be based on https://arxiv.org/pdf/1608.05854.pdf:
        #     Martin TB, Prunet S, Drissen L. Optimal fitting of Gaussian-apodized
        #     or under-resolved emission lines in Fourier transform spectra
        #     providing new insights on the velocity structure of NGC 6720. Monthly
        #     Notices of the Royal Astronomical Society. 2016 Sep 14;463(4):4223-38.
        # (thanks for the pointer Guillaume).
    end

    return krn
end

function smooth(vol::Array{<:Real,3},fwhm::Array{<:AbstractFloat,1})
    lim = ceil.(2*fwhm)
    x  = -lim[1]:lim[1]; x = smoothkern(fwhm[1],x); x  = reshape(x,(length(x),1,1))./sum(x)
    y  = -lim[2]:lim[2]; y = smoothkern(fwhm[2],y); y  = reshape(y,(1,length(y),1))./sum(y)
    z  = -lim[3]:lim[3]; z = smoothkern(fwhm[3],z); z  = reshape(z,(1,1,length(z)))./sum(z)

    vol1 = convsame(vol ,x)
    vol1 = convsame(vol1,y)
    vol1 = convsame(vol1,z)
    t    = typeof(vol)
    if t<:Array{<:Integer}
        vol1 = round.(vol1)
        vol1 = convert(t,vol1)
    end
end

"""

    A = spm_matrix(p)

Return an affine transformation matrix

 p[1]  - x translation
 p[2]  - y translation
 p[3]  - z translation
 p[4]  - x rotation about - {pitch} (radians)
 p[5]  - y rotation about - {roll}  (radians)
 p[6]  - z rotation about - {yaw}   (radians)
 p[7]  - x scaling
 p[8]  - y scaling
 p[9]  - z scaling
 p[10] - x affine
 p[11] - y affine
 p[12] - z affine

 A     - affine transformation matrix

---

 spm_matrix returns a matrix defining an orthogonal linear (translation,
 rotation, scaling or affine) transformation given a vector of
 parameters (p).  By default, the transformations are applied in the
 following order (i.e., the opposite to which they are specified):

 1) shear
 2) scale (zoom)
 3) rotation - yaw, roll & pitch
 4) translation

 SPM uses a PRE-multiplication format i.e. Y = A*X where X and Y are 4 x n
 matrices of n coordinates.
"""
function spm_matrix(p)

    #-Pad p with 'null' parameters
    #--------------------------------------------------------------------------
    q  = [0. 0. 0. 0. 0. 0. 1. 1. 1. 0. 0. 0.]
    p  = [p... q[(length(p) + 1):12]...]

    #-Translation / Rotation / Scale / Shear
    #--------------------------------------------------------------------------
    T  =   [1   0   0   p[1];
            0   1   0   p[2];
            0   0   1   p[3];
            0   0   0   1]

    R1  =  [1   0           0           0;
            0   cos(p[4])   sin(p[4])   0;
            0  -sin(p[4])   cos(p[4])   0;
            0   0           0           1]

    R2  =  [cos(p[5])   0   sin(p[5])   0;
            0           1   0           0;
           -sin(p[5])   0   cos(p[5])   0;
            0           0   0           1]

    R3  =  [cos(p[6])   sin(p[6])   0   0;
           -sin(p[6])   cos(p[6])   0   0;
            0           0           1   0;
            0           0           0   1]

    R   = R1*R2*R3

    Z   =  [p[7]   0       0       0;
            0      p[8]    0       0;
            0      0       p[9]    0;
            0      0       0       1]

    S   =  [1      p[10]   p[11]   0;
            0      1       p[12]   0;
            0      0       1       0;
            0      0       0       1]

    #-Affine transformation matrix
    #--------------------------------------------------------------------------
    A = T*R*Z*S
end

include("SPMlib.jl")

function objective(x::Array{<:Real}, G::Array{UInt8,3}, MG::Array{<:Real,2}, F::Array{UInt8,3}, MF::Array{<:Real,2},
                   s::Real=1., cf::String="nmi", fwhm::Array{<:Real,1}=[7., 7.])::AbstractFloat

    # Voxel sizes
    vxg  = sqrt.(sum(MG[1:3,1:3].^2,dims=1))
    sg   = s./vxg

    # Create the joint histogram
    H    = SPMlib.hist2(G,F, MF\spm_matrix(x)*MG ,sg)

    # Smooth the histogram
    lim  = ceil.(2.0.*fwhm);

    krn1 = smoothkern(fwhm[1],-lim[1]:lim[1]);
    krn1 = reshape(krn1,(length(krn1),1))./sum(krn1);
    H    = convn(H,krn1)

    krn2 = smoothkern(fwhm[2],-lim[2]:lim[2]);
    krn2 = reshape(krn2,(1,length(krn2)))./sum(krn2);
    H    = convn(H,krn2)

    # Compute cost function from histogram
    H  .+= eps(Float64)
    H  ./= sum(H)
    s1   = sum(H,dims=1)
    s2   = sum(H,dims=2)

    if cf=="mi"
        # Mutual Information:
        H  .*= log2.(H./(s2*s1))
        mi   = sum(H)
        o    = -mi
    elseif cf== "ecc"
        # Entropy Correlation Coefficient of:
        # Maes, Collignon, Vandermeulen, Marchal & Suetens (1997).
        # "Multimodality image registration by maximisation of mutual
        # information". IEEE Transactions on Medical Imaging 16(2):187-198
        H  .*= log2.(H./(s2*s1))
        mi   = sum(H)
        ecc  = -2*mi./(sum(s1.*log2.(s1)).+sum(s2.*log2.(s2)))
        o    = -ecc
    elseif cf=="nmi"
        # Normalised Mutual Information of:
        # Studholme,  Hill & Hawkes (1998).
        # "A normalized entropy measure of 3-D medical image alignment".
        # in Proc. Medical Imaging 1998, vol. 3338, San Diego, CA, pp. 132-143.
        nmi  = (sum(s1.*log2.(s1))+sum(s2.*log2.(s2)))/sum(H.*log2.(H))
        o    = -nmi
    elseif cf=="ncc"
        # Normalised Cross Correlation
        d    = size(H)
        i    = 1:d[1];
        j    = 1:d[2];
        m1   = sum(s2.*i')
        m2   = sum(s1.*j);
        sig1 = sqrt(sum(s2.*(i'.-m1).^2));
        sig2 = sqrt(sum(s1.*(j .-m2).^2));
        i    = (i.-m1).*ones(Float64,(1,d[2]))
        j    = (j.-m2)'.*ones(Float64,(d[1],1))
        ncc  = sum(sum(H.*i.*j))/(sig1*sig2)
        o    = -ncc
    else
        error("invalid cost function specified")
    end
    return o
end

include("Nii.jl")
include("Powell.jl")

function to8bit(nii)::Array{UInt8}
    mx    = maximum(nii.raw[:])
    mn    = minimum(nii.raw[:])
    if isa(nii.raw,Array{Int16}) || isa(nii.raw,Array{UInt16}) || isa(nii.raw,Array{UInt8}) || isa(nii.raw,Array{Int8})
        hlen = mx-mn+1
    else
        hlen = 2^15
    end
    si    = [1. mn; 1. mx]\[0.; hlen-1]
    f     = UInt16.(round.(nii.raw.*si[2].+si[1]))

    # Compute histogram
    h     = zeros(Int64,hlen);
    for i=1:length(f) h[f[i]+1]+=1; end

    # get lower and upper values from range that exclude a small number of outliers
    ch    = cumsum(h,dims=1)
    upper = findlast( ch.<=ch[end]*0.9999)-1
    lower = findfirst(ch.>=ch[end]*(1-0.9999))-1

    # Use these as min and max values for rescaling to UInt8
    mx    = (upper-si[1])/si[2]
    mn    = (lower-si[1])/si[2]
    si    = [1 0; 1 255]\[mn; mx]
    r     = (typeof(nii.raw)<:Array{<:Integer} ? rand(Float32,size(nii.raw)) : 0.0)
    f     = UInt8.(max.(min.(floor.((nii.raw.+r).*si[2].+si[1]),255f0),0f0))
end

function run(Pfixed::String,Pmoved::String, sep::Array{<:Real}=[4.0, 2.0, 1.0], cost::String="nmi", hsmo::Array{<:Real}=[5.,5.])
    seed!(100)

    Nfixed = Nii.nifti(Pfixed)
    Ffixed = Coreg.to8bit(Nfixed)
    Mfixed = Nii.vox2world(Nfixed)
    vx     = sqrt.(sum(Mfixed[1:3,1:3].^2,dims=1))[:]
#   fwhm   = sqrt.(max.(sep[end]^2 .- vx[:].^2, 0.0))./vx
#   Ffixed = smooth(Ffixed,fwhm)

    Nmoved = Nii.nifti(Pmoved)
    Fmoved = to8bit(Nmoved)
    Mmoved = Nii.vox2world(Nmoved)
    vx     = sqrt.(sum(Mmoved[1:3,1:3].^2,dims=1))[:]
#   fwhm   = sqrt.(max.(sep[end]^2 .- vx.^2, 0.0))./vx
#   Fmoved = smooth(Fmoved,fwhm)

    x     = zeros(6)
    o     = Inf
    tolsc = [0.02, 0.02, 0.02, 0.001, 0.001, 0.001]
    xi    = diagm(0=>20.0.*tolsc)

    for s in sep
       @printf("%s -> %s w. %s, %.3gmm samp. & %.3gx%.3g hist. smo.\n",
               splitdir(Pmoved)[end], splitdir(Pfixed)[end],
               cost, s, hsmo[1], hsmo[2])
        x,o = Powell.powell(x,xi,tolsc,objective,Ffixed, Mfixed, Fmoved, Mmoved, s, cost, hsmo)
    end
    return x,o
end

end # module

