module CertifiedHomotopyTracking

using Reexport
@reexport using Symbolics
@reexport using Nemo
@reexport using AbstractAlgebra
@reexport using LinearAlgebra
@reexport using MultivariatePolynomials
@reexport using DynamicPolynomials 
@reexport using GAP

export @setupfield

export straight_line_homotopy, specified_system, track

export track_complete_graph, get_permutations, str_convert, HCSystem, CompiledHomotopy, make_edge_system

# Source Code Include

# [Internals] -- helpers and utilities
include("internals/interval_arithmetic.jl")
include("internals/linear_algebra.jl")
include("internals/systems.jl")
include("internals/homotopy_constructor.jl") 
include("internals/taylor_model.jl") 
include("internals/moore_box.jl")       
include("internals/krawczyk.jl")        
include("internals/predictors.jl")
include("internals/tracking_modules.jl") 
include("internals/homogenize.jl") 
include("internals/complex_coefficients.jl")

# [Core] -- main functionalities
include("poly_setup.jl")
include("homotopy.jl")
include("results.jl")
include("tracking.jl")    
include("monodromy.jl")   

end # module
