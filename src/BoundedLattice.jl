export BoundedLattice, boundary, isperiodic

"""
    BoundedLattice{N, BC} <: AbstractLattice{N}

A lattice with boundary conditions `BC` on `N` dimension.
"""
abstract type BoundedLattice{N, BC <: BoundaryCondition} <: AbstractLattice{N} end


# Bounded Lattice Traits
"""
    boundary(lattice) -> BoundaryCondition

Returns lattice's boundary condition.
"""
function boundary end

"""
    isperiodic(lattice) -> Bool

Check if this lattice's boundary condition is periodic.
"""
function isperiodic end

boundary(::Type{<:BoundedLattice{N, BC}}) where {N, BC} = BC

isperiodic(::Type{<:BoundedLattice{N, Periodic}}) where N = true
isperiodic(::Type{<:BoundedLattice}) = false

# forward instance traits
isperiodic(::T) where {T <: BoundedLattice} = isperiodic(T)
boundary(::T) where {T <: BoundedLattice} = boundary(T)

nameof(::BoundedLattice) = "Bounded Lattice"

function show(io::IO, lattice::BoundedLattice)
    align = get(io, :align, 0)
    println(io, " " ^ align, nameof(lattice))
    println(io, " " ^ align, "boundary: ", boundary(lattice))
    println(io, " " ^ align, "size:     ", size(lattice))
end
