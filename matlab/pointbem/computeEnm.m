function Enm = computeEnm(b, epsIn, pqrData, q, Nmax)

Nq = length(q);

sphPoints = convertToSph(pqrData.xyz);
r = sphPoints(:,1);
phi = sphPoints(:,2);
theta = sphPoints(:,3);

Enm = zeros(Nmax,2*Nmax+1);

for n=0:Nmax-1
  iIndex=n+1;
  Pnm = legendre(n, cos(theta));

  for m=-n:n
	 jIndex=m+n+1;
	 ff = factorial(n-abs(m))/factorial(n+abs(m));
	 for k=1:Nq
		Enm(iIndex,jIndex) = Enm(iIndex,jIndex) + ff * q(k) * r(k)^n * Pnm(abs(m)+1,k) * exp(-i*m*phi(k));
	 end
  end
end
