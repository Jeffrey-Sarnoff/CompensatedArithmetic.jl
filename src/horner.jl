
#=
Accurate simple zeros of polynomials in floating point arithmetic
Stef Graillat
Computers and Mathematics with Applications 56 (2008) 1114–1120

   poly has N+1 coefficients a_0,a_1,..,a_N for x^0,x^1,..,x^N 
   
   poly = polynomial(coeffs[a_0,a_1,..,a_N])
   poly(x) = sum[i=0:1:N]( a_i * x^i )          this the the classical Horner scheme
   
   classicalHorner{R<:Real}(p::polynomial,x::R)
      n = length(p.coeffs)
      s = zeros(R, n)       # one-based indexing
      s[:] = p.coeffs[:]
      for i in (n-1):(-1):1
          s[i] += s[i+1]*x
      end
      s[1]
    end  

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

function sumHorner{T<:AbstractFloat}(a::Vector{T}, b::Vector{T}, x::T)
   vlen = min(length(a),length(b))
   r = a[vlen]+b[vlen]
   for i in (vlen-1): -1: 1
     r = r * x + (a[i]+b[i])
   end
   r
end

function csdHornerA{T<:AbstractFloat}(p::Poly{T}, x::T)
    r, pa, pb = eftHorner(p,x)
    c = sumHorner(pa,pb,x)
    r + c
end

# as accurate as if computed at twice the working precision and then rounded to working precision
function csdHorner{T<:AbstractFloat}(p::Poly{T}, x::T)
    r, pa, pb = eftHorner(p,x)
    pab = Poly(pa + pb)
    c = polyval(pab, x)
    r + c
end
