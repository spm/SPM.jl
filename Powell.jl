module Powell

using Printf

const gold  = (1.0+sqrt(5.0))/2.0 # Golden ratio
const gold1 = 1.0-(sqrt(5.0)-1.0)/2.0
const ϵ     = eps(Float64)
struct PF
    p::AbstractFloat
    f::AbstractFloat
    PF(p,f) = new(p,f)
    PF(a::Tuple{<:AbstractFloat,<:AbstractFloat}) = new(a[1],a[2])
end
struct FunType
    p::Array{<:AbstractFloat,1}
    p1::Array{<:AbstractFloat,1}
    func::Function
    args::Tuple
end


"""
Powell optimisation method

    FORMAT [p,f] = powell(p,xi,tolsc,func,varargin)
    p        - Starting parameter values
    xi       - columns containing directions in which to begin searching
    tolsc    - stopping criteria, optimisation stops when
                 sqrt(sum(((p-p_prev)./tolsc).^2))<1
    func     - name of evaluated function
    varargin - remaining arguments to func (after p)

    p        - final parameter estimates
    f        - function value at minimum

---

 Method is based on Powell's optimisation method described in
 Numerical Recipes (Press, Flannery, Teukolsky & Vetterling).

---
Copyright (C) 2001-2017 Wellcome Trust Centre for Neuroimaging
"""
function powell(p::Array{<:AbstractFloat,1},xi::Array{<:AbstractFloat,2},tolsc::Array{<:AbstractFloat},func,varargin...)
    f = func(p,varargin...)
    for iter=1:512
        #if numel(p)>1, @printf("iteration %d...\n", iter); end;            %-#
        ibig = length(p)
        pp   = p
        fp   = f
        del  = 0
        for i=1:length(p)
            ft = f
            p,junk,f = min1d(p,xi[:,i],func,f,tolsc,varargin)
            if abs(ft-f) > del
                del  = abs(ft-f)
                ibig = i
            end
        end
        if length(p)==1 || sqrt(sum(((p.-pp)./tolsc).^2))<1.0 || abs((f-fp)/(f+fp))<1e-6
            return p, f
        end
        ft = func(2.0.*p.-pp,varargin...)
        if ft < f
            p,xi[:,ibig],f = min1d(p,p.-pp,func,f,tolsc,varargin)
        end
    end
    return p, f
    #warning("Too many optimisation iterations")
end


"""
Line search for minimum
    (p,p1,f) = min1d(p,p1,func,f,tolsc,varargin)
"""
function min1d(p,p1,func,f,tolsc,varargin)
    global lnm = FunType(p,p1,func,varargin)
    tol    = 1.0/sqrt(sum((p1./tolsc).^2))
    t      = bracket(f)
    f,pmin = search(t,tol)
    p1     = p1.*pmin
    p      = p .+ p1

    if length(p)<12
        map(p->@printf("%-8.4g ",p),p)
        @printf("| %.5g\n", f)
    else
        @printf("%.5g\n", f)
    end

    return p, p1, f
end

"""
Reconstruct parameters and evaluate

    f = funeval(p)
"""
function funeval(p::AbstractFloat)
    global lnm # defined in min1d
    pt = lnm.p.+p.*lnm.p1
    f  = lnm.func(pt,lnm.args...)
   #println(p," ",f)
    return f
end

"""
Bracket the minimum (t[2]) between t[1] and t[3]
    t = bracket(f)
"""
function bracket(f)

    t = [PF(0.,f), PF(1.,funeval(1.)), PF(0.,f)]

    # if t[2] not better than t[1] then swap
    if t[2].f > t[1].f
        t[3] = t[1]
        t[1] = t[2]
        t[2] = t[3]
    end
    pp   = t[2].p + gold*(t[2].p-t[1].p)
    t[3] = PF(pp,funeval(pp))

    while t[2].f > t[3].f

        # fit a polynomial to t
        pol = ((map(e->e.p,t).-t[2].p).^(0.0:2.0)')\map(e->e.f,t)

        # minimum is when gradient of polynomial is zero
        # sign of pol(3) (the 2nd deriv) should be +ve
        if pol[3]>0
            # minimum is when gradient of polynomial is zero
            d    = -pol[2]/(2.0*pol[3]+ϵ)

            # A very conservative constraint on the displacement
            if d >  (1+gold)*(t[3].p-t[2].p)
                d = (1+gold)*(t[3].p-t[2].p)
            end
            pp  = t[2].p+d
        else
            # sign of pol[3] (the 2nd deriv) is not +ve
            # so extend out by golden ratio instead
            pp  = t[3].p+gold*(t[3].p-t[2].p)
        end

        # FUNCTION EVALUATION
        u = PF(pp,funeval(pp))

        if (t[2].p < u.p) == (u.p < t[3].p)

            # u is between t[2] and t[3]
            if u.f < t[3].f
                # minimum between t[2] and t[3] - done
                t[1] = t[2]
                t[2] = u
                return t
            elseif u.f > t[2].f
                # minimum between t[1] and u - done
                t[3] = u
                return t
            end
        end

        # Move all 3 points along
        t[1] = t[2]
        t[2] = t[3]
        t[3] = u
    end
    return t
end

"""
Brent's method for line searching - given that minimum is bracketed
    (f,p) = search(t, tol)
"""
function search(t, tol)

    # Current and previous displacements
    d     = Inf
    pd    = Inf

    # sort t into best first order
    t   = t[sortperm(map(e->e.f,t))];
    brk = [minimum(map(e->e.p,t)) maximum(map(e->e.p,t))];

    for iter=1:128
        # check stopping criterion
        if abs(t[1].p - 0.5*(brk[1]+brk[2]))+0.5*(brk[2]-brk[1]) <= 2.0*tol
            return t[1].f, t[1].p
        end

        # keep last two displacents
        ppd = pd
        pd  = d

        # fit a polynomial to t
        pol = ((map(e->e.p,t).-t[1].p).^(0.0:2.0)')\map(e->e.f,t)

        # minimum is when gradient of polynomial is zero
        d   = -pol[2]/(2*pol[3]+ϵ)
        pp  = t[1].p+d

        # check so that displacement is less than the last but two,
        # that the displaced point is between the brackets
        # and that the solution is a minimum rather than a maximum
        ϵ₂ = 2.0*ϵ*abs(t[1].p)+ϵ;
        if abs(d) > abs(ppd)/2 || pp < brk[1]+ϵ₂ || pp > brk[2]-ϵ₂ || pol[3]<=0
            # if criteria are not met, then golden search into the larger part
            d   = (t[1].p >= 0.5*(brk[1]+brk[2]) ? gold1*(brk[1]-t[1].p) : gold1*(brk[2]-t[1].p) )
            pp  = t[1].p+d
        end

        # FUNCTION EVALUATION
        u = PF(pp,funeval(pp))

        # Insert the new point into the appropriate position and update
        # the brackets if necessary
        #brk = (u.f<=t[1].f ? (u.p>=t[1].p ? [t[1].p, brk[2]] : [brk[1], t[1].p]) : (u.p < t[1].p ? [u.p, brk[2]] : [brk[1], u.p]))
        #t   = (u.f<=t[1].f ? [u,t[1],t[2]] : (u.f <= t[2].f ? [t[1],u,t[2]] : [t[1],t[2],u]))

        if u.f <= t[1].f
            if u.p >= t[1].p
                brk[1]=t[1].p
            else
                brk[2]=t[1].p
            end
            t[3] = t[2]
            t[2] = t[1]
            t[1] = u
        else
            if u.p < t[1].p
                brk[1]=u.p
            else
                brk[2]=u.p
            end
            if u.f <= t[2].f
                t[3] = t[2]
                t[2] = u
            elseif u.f <= t[3].f
                t[3] = u
            end
        end
    end
    return f,p
end

end # module
