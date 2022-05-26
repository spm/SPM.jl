module Hello
push!(LOAD_PATH, joinpath(pwd(),"src"))
import SPM

Base.@ccallable function julia_main(ARGS::Vector{String})::Cint

    if length(ARGS)<2
        println("Usage: coreg fixed.nii moved.nii samp cost hsmo1 hsmo2")
        println("  fixed - filename of the NIfTI image that remains fixed")
        println("  moved - filename of the NIfTI image that is moved")
        println("  samp  - sampling distance (default 2.0mm)")
        println("  cost  - cost function (default:\"nmi\", \"mi\", \"ecc\" or \"ncc\")")
        println("  hsmo1 - histogram smoothing (default:7.0)")
        println("  hsmo2 - histogram smoothing (default:7.0)")
        println("")
        return 1
    end

    samp    = [(length(ARGS)>=3 ? parse(Float64,ARGS[3]) : 2.)]
    cost    =  (length(ARGS)>=4 ? ARGS[4] : "nmi")
    hsmo    = [7., 7.]
    hsmo[1] =  (length(ARGS)>=5 ? parse(Float64,ARGS[5]) : 7.0)
    hsmo[2] =  (length(ARGS)>=6 ? parse(Float64,ARGS[6]) : hsmo[1])

    if length(ARGS)>6
        println("Too many input arguments.")
        return 1
    end

    try
        x,o = SPM.Coreg.run(ARGS[1],ARGS[2],samp,cost,hsmo)
        A   = SPM.Coreg.spm_matrix(x)
        println(A)
    catch
        println("An error occurred.")
        return 1
    end
    return 0
end

end

