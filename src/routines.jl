function insert_zero(x, loc)
    return x + x & ~(1 << (loc - 1) - 1)
end
