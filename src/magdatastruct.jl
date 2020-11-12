

##########################################

"""
$(TYPEDEF)

Structure containing the segments of a polygonal body.
To create an instance a set of indices have to be passed on.

# Fields 

$(TYPEDFIELDS)

"""
struct BodySegments2D
    "(x,y) for first set of vertices (beginning of segments)"
    ver1::SubArray 
    "(x,y) for second set of vertices (end of segments)"
    ver2::SubArray
    "total number of segments"
    nsegm::Integer

    function BodySegments2D(idx1::Union{Vector{<:Integer},UnitRange{<:Integer}},
                            vertices::Array{<:Real,2})
        # @assert ndims(ver1)==2
        @assert size(vertices,2)==2
        # circular shift to get second set of indices
        idx2 = circshift(idx1,-1)
        ## first set of vertices
        ver1 = view(vertices,idx1,:)
        ## second set of vertices
        ver2 = view(vertices,idx2,:)
        nsegm = size(ver1,1)
        return new(ver1,ver2,nsegm)
    end    
end

##########################################

"""
$(TYPEDEF)

Structure containing a set of polygonal bodies described by their segments and all vertices.
To create an instance, input an array of vectors of indices 
  (of vertices) for each body and the array of all the vertices.

# Fields 

$(TYPEDFIELDS)

"""
struct MagPolygBodies2D
    "array of bodies defined by their vertices"
    bo::Vector{BodySegments2D}
    "array of all vertices for all bodies"
    allvert::Array{<:Real,2}
 
    function MagPolygBodies2D(bodyindices::Vector{<:Vector{<:Integer}},allvert::Array{<:Real,2})
        @assert size(allvert,2)==2
        N=length(bodyindices)
        bodies = Vector{BodySegments2D}(undef,N)
        for i=1:N
            bodies[i] = BodySegments2D(bodyindices[i],allvert)
        end
        return new(bodies,allvert)
    end
end

##########################################

"""
$(TYPEDEF)

Structure containing the components of a magnetization vector, 
  i.e., module, inclination and declination angles.

# Fields 

$(TYPEDFIELDS)

"""
Base.@kwdef struct MagnetizVector
    "modulus "
    mod::Real
    "inclination in degrees"
    Ideg::Real
    "declination in degrees"
    Ddeg::Real
end

##########################################
