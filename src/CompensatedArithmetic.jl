module CompensatedArithmetic

!isdefined(Float) && typealias Float AbstractFloat;

export csdDot, csdDet2x2, csdCross3D, csdHorner

using ErrorfreeArithmetic


include("arith.jl")
include("det.jl")
include("horner.jl")
include("dotproduct.jl")

end # module
