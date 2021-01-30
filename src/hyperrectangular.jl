struct HyperCubic{N, B<:NTuple{N, AbstractBoundary}, T <: Real} <: AbstractLattice
    dims::NTuple{N, Int}
    bcs::B
    lc::T
    translation_vectors::Dict{Int, Vector{Coordinate{N, Int}}}

    function HyperCubic{N, B, T}(dims::NTuple{N, Int}, bcs::B, lc::T) where {N, B <: NTuple{N, AbstractBoundary}, T <: Real}
        check_boundaries(bcs)
        lc > 0 || throw(ArgumentError("Lattice constant must be positive!"))
        new{N, B, T}(dims, bcs, lc, Dict{Int, Vector{Coordinate{N, Int}}}())
    end
end
latticeconstants(l::HyperCubic{N}) where N = ntuple(_ -> l.lc, N)

struct HyperRectangular{N, B<:NTuple{N, AbstractBoundary}, T <: Real} <: AbstractLattice
    dims::NTuple{N, Int}
    bcs::B
    lcs::NTuple{N, T}  # lattice constants
    translation_vectors::Dict{Int, Vector{Coordinate{N, Int}}}

    function HyperRectangular{N, B, T}(dims::NTuple{N, Int}, bcs::B, lcs::NTuple{N, T}) where {N, B <: NTuple{N, AbstractBoundary}, T <: Real}
        check_boundaries(bcs)
        all(x -> x > 0, lcs) || throw(ArgumentError("All lattice constants must be positive!"))

        if all(x -> x == lcs[1], lcs)
            # hits the case where N == 1 as well, since we'd rather not have
            #  two different types of 1D lattices floating around...
            return HyperCubic{N, B, T}(dims, bcs, lcs[1])
        else
            return new{N, B, T}(dims, bcs, lcs, Dict{Int, Vector{Coordinate{N, Int}}}())
        end
    end
end

################
## Constructors
################

function HyperRectangular{N}(
        dims::Union{Int, NTuple{N, Int}},
        bcs::Union{AbstractBoundary, NTuple{N, AbstractBoundary}}=Periodic(),
        lcs::Union{T, NTuple{N, T}}=1.0
    ) where {N, T <: Real}
    # auto-expand singlet args into tuples
    dims = (dims isa NTuple{N, Int}) ? dims : ntuple(_ -> dims, N)
    bcs = (bcs isa NTuple{N, AbstractBoundary}) ? bcs : ntuple(_ -> bcs, N)
    lcs = (lcs isa NTuple{N, Real}) ? lcs : ntuple(_ -> lcs, N)
    return HyperRectangular{N, typeof(bcs), eltype(lcs)}(dims, bcs, lcs)
end

function HyperCubic{N}(
        dims::Union{Int, NTuple{N, Int}},
        bcs::Union{AbstractBoundary, NTuple{N, AbstractBoundary}}=Periodic(),
        lc::T=1.0
    ) where {N, T <: Real}
    # auto-expand singlet args into tuples
    dims = (dims isa NTuple{N, Int}) ? dims : ntuple(_ -> dims, N)
    bcs = (bcs isa NTuple{N, AbstractBoundary}) ? bcs : ntuple(_ -> bcs, N)
    return HyperCubic{N, typeof(bcs), T}(dims, bcs, lc)
end

for L in [:HyperCubic, :HyperRectangular]
    lcname = (L == :HyperCubic) ? :lc : :lcs

    # infer N from dims tuple
    @eval $L(dims::NTuple{N, Int}, bcs, $lcname) where N = $L{N}(dims, bcs, $lcname)
    @eval $L(dims::NTuple{N, Int}, bcs) where N = $L{N}(dims, bcs)
    @eval $L(dims::NTuple{N, Int}) where N = $L{N}(dims)

    @eval ndims(::$L{N}) where N = N

    @eval ==(a::$L{N, B}, b::$L{N, B}) where {N, B <: NTuple{N, AbstractBoundary}} = (
        dimensions(a) === dimensions(b)
        && boundaryconditions(a) === boundaryconditions(b)
        && latticeconstants(a) == latticeconstants(b)
        && a.translation_vectors == b.translation_vectors  # should we even be comparing the translation vectors?
    )
end

show(io::IO, l::HyperRectangular{N}) where N =
    println(io, "HyperRectangular{$N}($(l.dims), $(l.bcs), $(l.lcs))")

show(io::IO, l::HyperCubic{N}) where N =
    println(io, "HyperCubic{$N}($(l.dims), $(l.bcs), $(l.lc))")

const Chain = HyperCubic{1}
const Square = HyperCubic{2}
const Cubic = HyperCubic{3}

show(io::IO, l::Chain) = println(io, "Chain($(l.dims), $(l.bcs), $(l.lc))")
show(io::IO, l::Square) = println(io, "Square($(l.dims), $(l.bcs), $(l.lc))")
show(io::IO, l::Cubic) = println(io, "Cubic($(l.dims), $(l.bcs), $(l.lc))")

#############################
## Translation/Basis Vectors
#############################

primitivevectors(l::HyperRectangular{N, B, T}) where {N, B, T} = Diagonal{T}(collect(l.lcs))
metric(l::HyperRectangular{N, B, T}) where {N, B, T} = Diagonal{T}(collect(l.lcs) .^ 2)
ismetricdiag(::HyperRectangular) = true

primitivevectors(l::HyperCubic) = UniformScaling(l.lc)
metric(l::HyperCubic) = UniformScaling(l.lc ^ 2)
ismetricdiag(::HyperCubic) = true

function precompute_translation_vectors!(lattice::Chain, maximum_k::Int)
    d = lattice.translation_vectors
    for i in 1:maximum_k
        d[i] = [Coordinate(i), Coordinate(-i)]
    end
    return d
end

#############
## Neighbors
#############

# in 1D, Helical and Periodic BCs are equivalent
for BC in (:Periodic, :Helical), k in (:1, :K)
    @eval function neighbors(lattice::Chain{NTuple{1,$BC}}, site::Coordinate{1, Int}, ::Val{$k}) where K
        l = lattice.dims[1]
        c = site.coordinates[1]
        return [
            Coordinate(apply_boundary(Periodic(), l, c - $k)),
            Coordinate(apply_boundary(Periodic(), l, c + $k)),
        ]
    end
end

function neighbors(lattice::Chain{NTuple{1,Open}}, site::Coordinate{1, Int}, ::Val{K}) where K
    l = lattice.dims[1]
    sp, sm = site + Coordinate(K), site - Coordinate(K)
    if sp.coordinates[1] > l
        if sm.coordinates[1] < 1
            return []
        else
            return [sm]
        end
    else
        if sm.coordinates[1] < 1
            return [sp]
        else
            return [sm, sp]
        end
    end
end
