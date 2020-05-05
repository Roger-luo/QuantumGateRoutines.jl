using TupleTools

struct OrderedTuple{S <: Tuple}
    storage::S
end

function insert_zero(x, loc)
    return x + x & ~(1 << (loc - 1) - 1)
end

# assume locs is ordered
function insert_zeros(x, locs::NTuple)
    ptr = first(locs)
    prev = 1 << (ptr - 1) - 1
    ret = x & prev
    intervals = TupleTools.diff(locs)
    for (k, inc) in enumerate(intervals)
        ptr += inc - 1
        curr = 1 << (ptr - 1) - 1
        mask = curr & ~prev
        prev = curr
        ret += (x & mask) << k
    end
    ret += (x & ~prev) << length(locs)
    return ret
end

function instruct!(r::AbstractVector, ::Val{:X}, loc::Int)
    @inbounds for lhs in 0:(length(r) >> 1)-1
        # mask locations before
        # insert 0 at loc
        p = lhs + lhs & ~(1 << (loc - 1) - 1)
        # insert 1 at loc
        q = lhs + lhs & ~(1 << (loc - 1) - 1) + 1 << (loc - 1)
        # swap rows
        tmp = r[p + 1]
        r[p + 1] = r[q + 1]
        r[q + 1] = tmp
    end
    return r
end

function instruct!(r::AbstractVector, ::Val{:X}, locs::NTuple{N, Int}) where N
    mask = bmask(locs)
    fixed_loc = locs[end]
    @inbounds for lhs in 0:(length(r) >> 1) - 1
        p = insert_zero(lhs, fixed_loc)
        i = flip(p, mask) + 1
        j = p + 1

        tmp = r[i]
        r[j] = r[i]
        r[i] = tmp
    end
    return r
end

function instruct_x!(r::AbstractVector, loc::Int)
    pow2loc = 1 << (loc - 1)
    mask = ~(pow2loc - 1)
    @inbounds for ii in 0:(length(r) >> 1)-1
        # we will consider how to write this part later
        # insert 0 at loc
        p = ii + ii & mask
        # insert 1 at loc
        q = ii + ii & mask + pow2loc
        # swap rows
        tmp = r[p + 1]
        r[p + 1] = r[q + 1]
        r[q + 1] = tmp
    end
    return r
end

function threaded_instruct_x!(r::AbstractVector, loc::Int)
    Threads.@threads for ii in 0:(length(r) >> 1)-1
        # we will consider how to write this part later
        # insert 0 at loc
        p = ii + ii & ~(1 << (loc - 1) - 1)
        # insert 1 at loc
        q = ii + ii & ~(1 << (loc - 1) - 1) + 1 << (loc - 1)
        # swap rows
        @inbounds begin
            tmp = r[p + 1]
            r[p + 1] = r[q + 1]
            r[q + 1] = tmp
        end
    end
    return r
end

function benchmark(qubit_range)
    otimes = zeros(length(qubit_range)); threaded_times = zeros(length(qubit_range));
    for (i, k) in enumerate(qubit_range)
        r = rand(ComplexF64, 1<<k)
        t1 = @benchmark instruct_x!($r, 2)
        t2 = @benchmark threaded_instruct_x!($r, 2)
        otimes[i] = minimum(t1) |> time
        threaded_times[i] = minimum(t2) |> time 
    end
    return otimes, threaded_times
end

struct Pattern{P, N} end
Pattern(P::Symbol, N::Int) = Pattern{P, N}()
const patterns = [:diagonal, :anti_diagonal]

function routine_m(name::Symbol, ex::Expr)
    n = checksquare(ex)
    for each in patterns
        ex = transform(Pattern(each, n), name, ex)
        ex === nothing || return ex
    end
    return fallback(n, name, ex)
end

function checksquare(ex::Expr)
    ex.head === :vcat || throw(Meta.ParseError("Expect a matrix definition, got $ex"))
    m = length(ex.args)
    all(isequal(m), map(x->length(x.args), ex.args)) || throw(Meta.ParseError("Expect a square matrix, got $ex"))
    return m
end

function fallback(n, name, ex::Expr) end

transform(::Pattern, name, ex) = nothing
function transform(::Pattern{:diagonal, 2}, name, ex)
    M11, M12 = ex.args[1].args
    M21, M22 = ex.args[2].args
    iszero(M21) && iszero(M12) || return

    quote
        function instruct!(r::AbstractVector, ::Val{$(QuoteNode(name))}, locs::NTuple{N, Int})
            mask = bmask(locs)
            n = length(locs)
            for b in 0:length(r)-1
                m = count_ones(b & mask)
                α = literal_pow($(Val(M11)), n-m) * literal_pow($(Val(M22)), m)
                mulrow!(r, b+1, α)
            end
        end
    end
end
