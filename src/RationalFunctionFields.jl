module RationalFunctionFields

# General packages
using Combinatorics
using LinearAlgebra
using Random
using TimerOutputs

# Algebra packages
using AbstractAlgebra
using Nemo
using Groebner
using ParamPunPam
using ParamPunPam: reduce_mod_p!, specialize_mod_p, AbstractBlackboxIdeal
ParamPunPam.enable_progressbar(false)


export RationalFunctionField

# helper functions
export generators

# membership functions
export check_algebraicity, field_contains, issubfield, fields_equal
export constructive_membership

# simplification
export simplified_generating_set

export IdealMQS, fractionfree_generators_raw

const _to = TimerOutputs.TimerOutput()

include("util.jl")
include("IdealMQS.jl")
include("Field.jl")
include("rankings.jl")
include("normalforms.jl")
include("simplification.jl")
include("constructive_membership.jl")

end
