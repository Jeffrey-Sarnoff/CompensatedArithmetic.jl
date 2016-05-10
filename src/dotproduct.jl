#=
    ref: Ogita, Rump, Oishi (SISC 05)
=#

""" dot(x,y) as if computed at twice working precision """
function csdDot{T<:AbstractFloat}(x::Vector{T}, y::Vector{T})
   s=c=p=q=r=zero(T)
   n = length(x)
   if ((n == 0)|(n != length(y)))
      throw(ErrorException("vectors must be of the same nonzero length"))
   end
   s,c = eftMul(x[1], y[1])
   for i in 2:n
       p,q = eftMul(x[i], y[i])
       s,r = eftAdd(s, p)
       c += q+r
   end
   s+c
end   
