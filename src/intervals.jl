"""
    within_interval(x, interval)

Return the indices of the elements of `x` that lie within the given
`interval`.
"""
within_interval(x::AbstractRange, interval::Interval) = findall(in(interval), x)
