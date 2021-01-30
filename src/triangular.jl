struct Triangular{B <: NTuple{2, AbstractBoundary}} <: AbstractLattice
    dims::NTuple{2, Int}
    bcs::B
    translation_vectors::Dict{Int, Vector{Coordinate{2, Int}}}
    function Triangular{B}(dims::NTuple{2, Int}, bcs::B) where B <: NTuple{2, AbstractBoundary}
        check_boundaries(bcs)
        new{B}(dims, bcs, Dict{Int, Vector{Coordinate{2, Int}}}())
    end
end


for L in [:Honeycomb, :Kagome]

    @eval struct $L{B <: NTuple{2, AbstractBoundary}} <: AbstractLattice
        dims::NTuple{2, Int}
        bcs::B
        translation_vectors::Dict{Tuple{Int, Int}, Vector{Coordinate{3, Int}}}
        function $L{B}(dims::NTuple{2, Int}, bcs::B) where B <: NTuple{2, AbstractBoundary}
            check_boundaries(bcs)
            new{B}(dims, bcs, Dict{Tuple{Int, Int}, Vector{Coordinate{3, Int}}}())
        end
    end

    @eval hasunitcell(::$L) = Val{true}()
end


for L in [:Triangular, :Honeycomb, :Kagome]

    @eval ==(a::$L{B}, b::$L{B}) where B <: NTuple{2, AbstractBoundary} = (
        dimensions(a) === dimensions(b)
        && boundaryconditions(a) === boundaryconditions(b)
        && a.translation_vectors == b.translation_vectors
    )

    @eval ndims(::$L) = 2

    @eval primitivevectors(::$L) = [1.0 cosd(60); 0.0 sind(60)]
    @eval metric(::$L) = Symmetric([1.0 0.5; 0.5 1.0])
    @eval ismetricdiag(::$L) = false

    @eval $L(dims::NTuple{2, Int}, bcs::NTuple{2, AbstractBoundary}) = $L{typeof(bcs)}(dims, bcs)
    @eval $L(dims::NTuple{2, Int}, bc::AbstractBoundary=Periodic()) = $L(dims, ntuple(x->bc, 2))

end


#############################
## Unitcell Properties
#############################

unitcellsize(::Honeycomb) = 2
unitcellsize(::Kagome) = 3

unitcell_translation_vectors(::Honeycomb) = [
    0.0 1/3;
    0.0 1/3
]

unitcell_translation_vectors(::Kagome) = [
    0.0 0.5 0.0;
    0.0 0.0 0.5
]
