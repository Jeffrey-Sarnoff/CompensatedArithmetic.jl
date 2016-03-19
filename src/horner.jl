# shortened to return 2 values instead of three
function eftHorner(p,x)
   deg = degree(p); n=deg+1
   zHi=0.0;zLo=0.0
   r = coeffs(p)[n]
   for i in deg:(-1):1
      zHi,zLo = eftMul(r,x)
      r = zHi + coeffs(p)[i]
    end
    zHi,zLo
end

