"""
    within_interval(x, interval)

Return the indices of the elements of `x` that lie within the given
`interval`.
"""
function within_interval(x::AbstractRange, interval::Interval{L,R}) where {L,R}
    reversed = if step(x) < 0
        x = reverse(x)
        true
    else
        false
    end
    N = length(x)
    δx = step(x)

    l = leftendpoint(interval)
    r = rightendpoint(interval)

    x₁ = x[1]

    a = max(ceil(Int, (l-x₁)/δx) + 1, 1)
    a += (a ≤ N && (x[a] == l && L == :open || x[a] < l))

    b = min(ceil(Int, (r-x₁)/δx) + 1, N)
    b -= (b ≥ 1 && (x[b] == r && R == :open || x[b] > r))

    if reversed
        a,b = b,a
        a = N - a + 1
        b = N - b + 1
    end

    a:b
end
