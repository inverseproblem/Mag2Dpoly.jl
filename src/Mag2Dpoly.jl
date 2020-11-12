"""
Mag2Dpoly

A module to perform magnetic anomaly calculations for 2D polygonal bodies.

# Exports

$(EXPORTS)
"""
module Mag2Dpoly

using DocStringExtensions

export BodySegments2D,MagPolygBodies2D,MagnetizVector
export tmagpolybodies2D,tmagpolybodies2Dgen
export tmagpoly2D,tmagpoly2Dgen

include("magutils.jl")
include("magdatastruct.jl")
include("mag2dpolybodies.jl")


end # module
