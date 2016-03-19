# William Kahan's method, accuracy is within 2 ulp

function csd_ad_minus_bc{T<:AbstractFloat}(a::T, b::T, c::T, d::T)  
    bcHi  = b*c
    bcLo  = fma(-b, c,  bcHi)  # exact
    hi    = fma( a, d, -bcHi)  # close
    hi - bcLo                  # compensated
end

csd_ab_minus_cd{T<:AbstractFloat}(a::T, b::T, c::T, d::T) = ad_minus_bc(a,c,d,b)

function csdDet2x2{T<:AbstractFloat}(m2x2::Matrix{T})
   if size(k) != (2,2)
       throw(DomainError())
   end
   csd_ad_minus_bc(reshape(m2x2,4,1)...)
end

function csdCross3D{T<:AbstractFloat}(a::Vector{T}, b::Vector{T})
    a1,a2,a3 = a[1],a[2],a[3]
    b1,b2,b3 = b[1],b[2],b[3]
    
    x = ad_minus_bc(a2, a3, b2, b3)
    y = ad_minus_bc(a3, a1, b3, b1)
    z = ad_minus_bc(a1, a2, b1, b2)
    
    [x,y,z]
end
