
#=
Algorithms for Accurate, Validated and Fast Polynomial Evaluation
Stef Graillat, Philippe Langlois and Nicolas Louvet
Japan J. Indust. Appl. Math., 26 (2009), 191–214

eftHorner(p_n,x) computes
(i)  the floating point evaluation Horner(p_n,x)
(ii) two polynomials pa, pb, of degree n-1
such that
  Horner(p_n,x), pa, pb = eftHorner(p_n,x)
  p_n(x) == Horner(p_n,x) + (pa + pb)(x)
so
  eftHorner is an errorfree transformation for polynomial evaluation with the Horner algorithm
     (unless there is underflow)
=#

function eftHorner{T<:AbstractFloat}(p::Poly{T},x::T)
   deg = degree(p); n=deg+1
   cofs = coeffs(p)
   pa = zeros(T,n)
   pb = zeros(T,n)
   r = cofs[n]
   for i in deg:(-1):1
      t, pa[i] = eftMul(r, x)
      r, pb[i] = eftAdd(t, cofs[i])
    end
    r, pa, pb
end

# as accurate as if computed at twice the working precision and then rounded to working precision
function csdHorner{T<:AbstractFloat}(p::Poly{T}, x::T)
    r, pa, pb = eftHorner(p,x)
    pab = Poly(pa + pb)
    c = polyval(pab, x)
    r + c
end
