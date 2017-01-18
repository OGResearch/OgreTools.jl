__precompile__()

module OgreTools

depsjl = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
isfile(depsjl) ? include(depsjl) : error("OgreTools not properly ",
    "installed. Please run\nPkg.build(\"OgreTool\")")

include("ModelParser.jl")


end # module
