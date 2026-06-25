module CertifiedHomotopyTracking

using Reexport
@reexport using Symbolics
@reexport using Nemo
@reexport using AbstractAlgebra
@reexport using LinearAlgebra
@reexport using MultivariatePolynomials
@reexport using DynamicPolynomials 
@reexport using GAP

export straight_line_homotopy

export get_permutations, SpecializedHomotopy,
       CompiledHomotopy, HomotopySourceData, make_edge_system, MonodromyResult
export collect_hc_trace
export compile_system, AlgebraicVarietySystem, VarietyBox,
       VarietyApproximation,
       variety_system, system, evaluate_system, jacobian_system,
       certified_variety_approximation, export_variety_obj
export PathBox, PathVisualization, path_boxes, export_path_tikz, export_path_obj

# Source Code Include

# [Internals] -- helpers and utilities
include("internals/interval_arithmetic.jl")
include("internals/linear_algebra.jl")
include("internals/complex_coefficients.jl")
include("internals/systems.jl")
include("internals/homotopy_constructor.jl") 
include("internals/taylor_model.jl") 
include("internals/krawczyk.jl")        
include("threading.jl")
include("variety_system.jl")
include("internals/moore_box.jl")       
include("internals/predictors.jl")
include("internals/tracking_modules.jl") 
include("internals/homogenize.jl") 

# [Core] -- main functionalities
include("poly_setup.jl")
include("homotopy.jl")
include("visualization.jl")
include("results.jl")
include("tracking.jl")    
include("monodromy.jl")   
include("hc_trace.jl")
include("posteriori_certification.jl")
include("posteriori_tracker.jl")

end # module
