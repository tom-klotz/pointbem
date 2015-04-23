function bderiv = besselDerivative(bd0, M, z, b)
bderiv = zeros(M+1,1);
bderiv(1) = bd0;
for n=1:M
    bderiv(n+1) = (1/(2*n+1))*(n*b(1+n-1)-(n+1)*b(1+n+1));
end
