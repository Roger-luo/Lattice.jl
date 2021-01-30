"""
    AbstractLattice

Abstract type for general lattices.
For a more concrete definition please refer the following material:

- Lattice (group): https://en.wikipedia.org/wiki/Lattice_(group)
- Lattice Model: https://en.wikipedia.org/wiki/Lattice_model_(physics)
- Lattice Graph: https://en.wikipedia.org/wiki/Lattice_graph
"""
abstract type AbstractLattice end

unitcellsize(l::AbstractLattice) = 1
hasunitcell(::AbstractLattice) = Val{false}()

_ncoords(::Val{true}, l::AbstractLattice) = 1 + ndims(l)
_ncoords(::Val{false}, l::AbstractLattice) = ndims(l)
ncoordinates(l::AbstractLattice) = _ncoords(hasunitcell(l), l)

latticeconstants(l::AbstractLattice) = l.lcs
boundaryconditions(l::AbstractLattice) = l.bcs
dimensions(l::AbstractLattice) = l.dims
size(l::AbstractLattice) = dimensions(l)
size(l::AbstractLattice, d::Int) = size(l)[d]

primitivevectors(::AbstractLattice) = error("Not implemented")
metric(l::AbstractLattice) = (A = primitivevectors(l); A' * A)
ismetricdiag(l::AbstractLattice) = isdiag(metric(l))

==(a::L, b::L) where {L <: AbstractLattice} = false





function size(::Val{true}, l::AbstractLattice, dim::Int)
    if dim > ncoordinates(l)
        return 1  # copying behaviour of `size` for multidimensional arrays
    elseif dim == ncoordinates(l)
        return unitcellsize(l)
    else
        return l.dims[dim]
    end
end
size(::Val{false}, l::AbstractLattice, dim::Int) = (dim > ndims(l)) ? 1 : l.dims[dim]
size(l::AbstractLattice, dim::Int) = size(hasunitcell(l), l, dim)

size(::Val{true}, l::AbstractLattice) = (l.dims..., unitcellsize(l))
size(::Val{false}, l::AbstractLattice) = l.dims
size(l::AbstractLattice) = size(hasunitcell(l), l)

axes(l::AbstractLattice, dim::Int) = Base.OneTo(size(l, dim))
axes(l::AbstractLattice) = map(d -> Base.OneTo(d), size(l))

strides(l::AbstractLattice) = (1, cumprod(size(l)[1:end-1])...)
nsites(l::AbstractLattice) = prod(size(l))


struct WeightedLattice{L <: AbstractLattice, W} <: AbstractLattice
    lattice::L
    weights::W
end

struct CoordinateLattice{L <: AbstractLattice, Coordinate} <: AbstractLattice
    lattice::L
    coordinates::Vector{Coordinate}
end

struct GraphLattice{Graph} <: AbstractLattice
    graph::Graph
end


# struct HasCoordinate end

# function coordinate(::HoneyComb, id)
# end

# function coordinate(lattice, id)
# end

# function coordinate(::HasCoordinate, lattice, id)
# end

# for (i, j) in edges(lattice; bond=k)
# end

# # <=>
# for i in sites(lattice)
#     for j in neighbors(lattice, i, k)
#         if i < j
#             value = get_term(H, i, j)
#         end
#     end
# end

# struct SitesIt{L}
#     lattice::L
# end

# struct BondsIt{L}
#     lattice::L
#     bond::Int
# end

# function sites end
# function neighbors end
# function bonds end

# Base.iterate(it::BondsIt, st) = iterate_lattice(it, it.lattice, st)
# function iterate_lattice(::BondsIt, lattice, st) end
