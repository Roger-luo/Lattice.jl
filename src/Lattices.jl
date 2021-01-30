module Lattices

import Base: size, ndims, length, show, nameof, ==, axes, strides, checkbounds
using Base.Iterators
using LinearAlgebra

export AbstractLattice, WeightedLattice, CoordinateLattice

export neighbors, bonds, sites, translation_vectors, to_site_id, to_coordinate, nsites

export HyperRectangular, HyperCubic, Chain, Square, Cubic
export Triangular, Honeycomb, Kagome

export AbstractBoundary, Periodic, Open, Helical

export Coordinate

include("coordinate.jl")
include("boundaries.jl")
include("types.jl")
include("neighbors.jl")
include("sites.jl")
include("bonds.jl")
include("hyperrectangular.jl")
include("triangular.jl")


#######################
## Boundary Conditions
#######################

function _apply_boundary_conditions(::Val{true}, bcs::NTuple{M, Periodic}, dims::NTuple{M, Int}, site::Coordinate{N, Int})::Coordinate{N, Int} where {M, N}
    return Coordinate(apply_boundary.(bcs, dims, site.coordinates[1:end-1])..., site.coordinates[end])
end

function _apply_boundary_conditions(::Val{true}, bcs::NTuple{M, B}, dims::NTuple{M, Int}, site::Coordinate{N, Int})::Union{Coordinate{N, Int}, Nothing} where {M, N, B <: AbstractBoundary}
    coords = site.coordinates
    clipped = _apply_boundary_conditions(Val{false}(), bcs, dims, Coordinate(coords[1:end-1]))
    if clipped === nothing
        return nothing
    else
        return Coordinate(clipped.coordinates..., coords[end])
    end
end

function _apply_boundary_conditions(::Val{false}, bcs::NTuple{N, Periodic}, dims::NTuple{N, Int}, site::Coordinate{N, Int})::Coordinate{N, Int} where N
    return Coordinate(apply_boundary.(bcs, dims, site.coordinates)...)
end

function _apply_boundary_conditions(::Val{false}, bcs::NTuple{N, B}, dims::NTuple{N, Int}, site::Coordinate{N, Int})::Union{Coordinate{N, Int}, Nothing} where {N, B <: AbstractBoundary}
    coords = site.coordinates
    clipped = zeros(Int, N)
    for i in 1:N
        c = apply_boundary(bcs[i], dims[i], coords[i])
        if isnothing(c)
            return nothing
        else
            clipped[i] = c
        end
    end
    return Coordinate(clipped...)
end



for L in [:Triangular, :Honeycomb, :Kagome, :HyperCubic]

    @eval function apply_boundary_conditions(lattice::$L, site::Coordinate)
        return _apply_boundary_conditions(hasunitcell(lattice), lattice.bcs, lattice.dims, site)
    end

end


end
