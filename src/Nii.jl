"""
A wrapper for the NIfTI module
"""
module Nii

export nifti,vox2world

using NIfTI
using LinearAlgebra

const ni_units = Dict([
     ( 0,(   1,"UNKNOWN")),
     ( 1,(1000,"m")),
     ( 2,(   1,"mm")),
     ( 3,(1e-3,"um")),
     ( 8,(   1,"s")),
     (16,(1e-3,"ms")),
     (24,(1e-6,"us")),
     (32,(   1,"Hz")),
     (40,(   1,"ppm")),
     (48,(   1,"rads"))]);

# Convert between first voxel at [1,1,1] (Julia) and first voxel at [0,0,0] (NIfTI header)
const One_off = [1.0 0.0 0.0 -1.0; 0.0 1.0 0.0 -1.0; 0.0 0.0 1.0 -1.0; 0.0 0.0 0.0 1.0]

function nifti(fname)
    return(niread(fname))
end

"""
    M = q2m(q)

Generate a rotation matrix from a quaternion xi+yj+zk+w,
where q = [x y z], and w = 1-x^2-y^2-z^2.
See: http://skal.planet-d.net/demo/matrixfaq.htm
"""
function q2m(q)
    q = Float64.(q[1:3]) # Assume rigid body
    w = sqrt(1 - sum(q.^2))
    x = q[1]
    y = q[2]
    z = q[3]
    if w<1e-7
        w = 1/sqrt(x*x+y*y+z*z)
        x = x*w
        y = y*w
        z = z*w
        w = 0
    end
    xx = x*x; yy = y*y; zz = z*z; ww = w*w;
    xy = x*y; xz = x*z; xw = x*w;
    yz = y*z; yw = y*w; zw = z*w;
    M = [ (xx-yy-zz+ww)      2*(xy-zw)      2*(xz+yw) 0
              2*(xy+zw) (-xx+yy-zz+ww)      2*(yz-xw) 0
              2*(xz-yw)      2*(yz+xw) (-xx-yy+zz+ww) 0
                     0              0              0  1]
end

"""
    q = m2q(M)

Convert from rotation matrix to quaternion form
See: http://skal.planet-d.net/demo/matrixfaq.htm
"""
function m2q(M)
    M = Float64.(M)
    d = diag(M[1:3,1:3])
    t = sum(d) + 1
    if t>0.5
        s = sqrt(t)*2;
        q = [(M[3,2]-M[2,3])/s, (M[1,3]-M[3,1])/s, (M[2,1]-M[1,2])/s, 0.25*s];
    else
        t = findmax(d);
        t = t[2];
        if t==1
            s = 2*sqrt(1 + M[1,1] - M[2,2] - M[3,3]);
            q = [0.25*s (M[1,2]+M[2,1])/s (M[3,1]+M[1,3])/s (M[3,2]-M[2,3])/s]';
        elseif t==2
            s = 2*sqrt(1 + M[2,2] - M[1,1] - M[3,3]);
            q = [(M[1,2]+M[2,1])/s 0.25*s (M[2,3]+M[3,2])/s (M[1,3]-M[3,1])/s ]';
        else
            s = 2*sqrt(1 + M[3,3] - M[1,1] - M[2,2]);
            q = [(M[3,1]+M[1,3])/s (M[2,3]+M[3,2])/s 0.25*s (M[2,1]-M[1,2])/s]';
        end
    end
    if q[4]<0 q = -q end # w must be +ve
    q = Float32.(q)
end


"""
Get the qform matrix from the NIfTI header (0-indexed)
"""
function qmat(hdr::NIfTI.NIfTI1Header)

    I₄₃ = Matrix{Float64}(I,4,3)

    # Rotations from quaternions
    R = q2m([hdr.quatern_b, hdr.quatern_c, hdr.quatern_d])

    # Translations
    T = [I₄₃ Float64.([hdr.qoffset_x, hdr.qoffset_y, hdr.qoffset_z, 1])]

    # Zooms.  Note that flips are derived from the first
    # element of pixdim, which is normally unused.
    n = min(hdr.dim[1],3);
    Z = [hdr.pixdim[2:(n+1)]..., ones(Float64,4-n)...];
    Z[Z.<0.0] .= 1.0;
    if hdr.pixdim[1]<0 Z[3] = -Z[3] end
    Z = diagm(0 => Z)
    return(T*R*Z)
end


"""
Get the sform matrix from the NIfTI header (0-indexed)
"""
function smat(hdr::NIfTI.NIfTI1Header)
    # read s-form matrix
    M  = [[hdr.srow_x...]'; [hdr.srow_y...]'; [hdr.srow_z...]'; 0.0 0.0 0.0 1.0];
end


"""
Get the voxel-to-world mapping matrix from the NIfTI header (1-indexed)
"""
function vox2world(hdr::NIfTI.NIfTI1Header)

    # Convert between first voxel at [1,1,1] (Julia) and first voxel at [0,0,0] (NIfTI header)
    One_off = [1.0 0.0 0.0 -1.0; 0.0 1.0 0.0 -1.0; 0.0 0.0 1.0 -1.0; 0.0 0.0 0.0 1.0]

    # Convert units to mm
    s,units = get(ni_units, hdr.xyzt_units & 7, (1, "UNKNOWN"));
    U       = diagm(0 => [s,s,s,1.0])

    if hdr.sform_code!=0
        M = U*smat(hdr)*One_off
    elseif hdr.qform_code!=0
        M = U*qmat(hdr)*One_off
    else
        error("missing positional information")
    end
    return(M)
end


"""
Get the voxel-to-world mapping matrix from the NIfTI object (1-indexed)
"""
function vox2world(nii::NIVolume)
    vox2world(nii.header)
end

end # module

