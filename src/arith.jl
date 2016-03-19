function csdDot{T<:AbstractFloat}(x::Vector{T}, y::Vector{T})
    xlen = length(x)
    if xlen != length(y)
       throw(ErrorException("vector lengths differ")
    elseif xlen==0
       return(x)
    end
    
    sHi, compensation = eftMul(x[1], y[1])
    for i in 2:xlen
        pHi, pLo = eftMul(x[i], y[i])
        sHi, sLo = eftAdd(sHi, pHi)
        compensation += pLo + sLo
    end
    
    sHi+compensation                          # as if rounded from 2x significant bits
end