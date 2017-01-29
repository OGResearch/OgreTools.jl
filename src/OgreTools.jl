#__precompile__()

module OgreTools

# check if the package was properly built
depsjl = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
isfile(depsjl) ? include(depsjl) : error("OgreTools not properly ",
    "installed. Please run\nPkg.build(\"OgreTools\")")

# import Base methods to be overloaded
import Base: show, copy, print, ==

# include definition and methods of the ParsedModel type
include("ParsedModel.jl")
export ParsedModel, parseFile

# include definition and methods of the Date type
include("Date.jl")
export Date, mm, qq, yy, ypf, dat2str

# include definition and methods of the TimeSeries type
include("TimeSeries.jl")
export TimeSeries,
       pct, apct, hpf,
       extend!, justify!

# include definition and methods of the DataBase type
include("DataBase.jl")
export DataBase, justify!, db2array,
       calc_irf

# include definition and methods of the SolveModel type
include("SolveModel.jl")
export SolveModel, solve!,
       get_endog_ind, get_stacked_jac,
       eval_equations, eval_derivatives,
       residual, residual!, jacobian,
       dac!, homotopy!,
       put_x_to_db!, get_x_from_db,
       endogenize!, exogenize!,
       shock!, get_result, check_model,
       check_resid

end # module
