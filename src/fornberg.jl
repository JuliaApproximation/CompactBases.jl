function fornberg_weights!(c::AbstractMatrix{T}, x, ξ) where T
    nd,m = size(c) .- 1
    c = OffsetArray(c, 0:nd, 0:m)
    x = OffsetArray(x, 0:nd)
    length(x) == nd+1 ||
        throw(DimensionMismatch("Number of grid points does not match requested number of derivative coefficients"))

    c₁ = one(T)
    c₄ = x[0] - ξ

    c .= zero(T)
    c[0,0] = 1

    for i = 1:nd
        mn = min(i,m)
        c₂ = one(T)
        c₅ = c₄
        c₄ = x[i] - ξ
        for j = 0:i-1
            c₃ = x[i] - x[j]
            c₂ *= c₃
            for k = mn:-1:1
                c[i,k] = c₁*(k*c[i-1,k-1] - c₅*c[i-1,k])/c₂
            end
            c[i,0] = -c₁*c₅*c[i-1,0]/c₂
            for k = mn:-1:1
                c[j,k] = (c₄*c[j,k] - k*c[j,k-1])/c₃
            end
            c[j,0] = c₄*c[j,0]/c₃
        end
        c₁ = c₂
    end
    c
end

function fornberg_weights(x, ξ::T, m::Integer) where T
    c = zeros(promote_type(eltype(x),T), length(x), m+1)
    fornberg_weights!(c, x, ξ)
    c
end

function fornberg_all_weights!(c::AbstractArray{T,3}, x, ξ) where T
    nd,m = size(c,1) .- 1,size(c,3) .- 1
    c = OffsetArray(c, 0:nd, 0:nd, 0:m)
    x = OffsetArray(x, 0:nd)
    length(x) == nd+1 ||
        throw(DimensionMismatch("Number of grid points does not match requested number of derivative coefficients"))

    c .= zero(T)
    c[0,0,0] = 1

    c₁ = one(T)
    c₄ = x[0] - ξ

    for i = 1:nd
        mn = min(i,m)
        c₂ = one(T)
        c₅ = c₄
        c₄ = x[i] - ξ
        for j = 0:i-1
            c₃ = x[i] - x[j]
            c₂ *= c₃
            c[j,i,0] = c₄*c[j,i-1,0]/c₃
            for k = 1:mn
                c[j,i,k] = (c₄*c[j,i-1,k] - k*c[j,i-1,k-1])/c₃
            end
        end
        c[i,i,0] = -c₁*c₅*c[i-1,i-1,0]/c₂

        for k = 1:mn
            c[i,i,k] = c₁*(k*c[i-1,i-1,k-1] - c₅*c[i-1,i-1,k])/c₂
        end
        c₁ = c₂
    end
    c
end

function fornberg_all_weights(x, ξ::T, m::Integer) where T
    c = zeros(promote_type(eltype(x),T), length(x), length(x), m+1)
    fornberg_all_weights!(c, x, ξ)
    c
end

function fornberg_test_driver(x, ξ)
    pnum(n::AbstractFloat) = printfmt("{1:12.4f}", n)
    pnum(n) = print(format(n, width=12))
    pnum(n::Rational) =
        print(format(isinteger(n) ? Int(n) : n, width=12))

    M = min(4,length(x)-1)
    C = fornberg_all_weights(x, ξ, M)
    N = size(C,2)-1
    print("j ")
    for j = 0:N
        pnum(j)
    end
    println()
    print("x ")
    foreach(pnum, x)
    println()
    println()
    for (ki,k) = enumerate(0:M)
        printfmtln("k = {1:2d}", k)
        for i = k:N
            ii = i+1
            printfmt("{1:2d}", i)
            foreach(pnum, view(C,1:ii,ii,ki))
            println()
        end
        println()
    end
    println()
    C = fornberg_weights(x, ξ, M)
    for (ki,k) = enumerate(0:M)
        printfmt("{1:2d}", k)
        foreach(pnum, view(C,:,ki))
        println(stdout)
    end
end
