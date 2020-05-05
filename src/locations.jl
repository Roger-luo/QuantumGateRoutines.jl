using TupleTools

abstract type AbstractLocations end

"""
    SortedLocations{S} <: AbstractLocations

Type for sorted location, where `S` is the storage type, it
can be any iteratable returns an integer.
"""
struct SortedLocations{S} <: AbstractLocations
    storage::S
end

"""
    RandomLocations{S} <: AbstractLocations

Type for random locations, where `S` is the storage type, it can
be any iteratable return an integer.
"""
struct RandomLocations{S} <: AbstractLocations
    storage::S
end

# forward
Base.length(x::AbstractLocations) = length(x.storage)
Base.getindex(x::AbstractLocations, k::Int) = x.storage[k]
Base.iterate(x::AbstractLocations) = iterate(x.storage)
Base.iterate(x::AbstractLocations, st) = iterate(x.storage, st)
Base.eltype(x::AbstractLocations) = eltype(x.storage)

Base.sort(x::SortedLocations; kwargs...) = x
Base.sort(x::RandomLocations; kwargs...) = SortedLocations(sort(x.storage; kwargs...))
Base.sort(x::RandomLocations{<:Tuple}; kwargs...) = SortedLocations(TupleTools.sort(x.storage); kwargs...)
