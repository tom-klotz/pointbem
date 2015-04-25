function phiReac = computePot(Bnm, pqrData, Nmax)

Nq = length(pqrData.q);
phiReac = zeros(Nq,1);

sphPoints = convertToSph(pqrData.xyz);
r     = sphPoints(:,1);
phi   = sphPoints(:,2);
theta = sphPoints(:,3);

for n=0:Nmax-1
  iIndex=n+1;
  Pnm = legendre(n, cos(theta));
  
  for m=-n:n
	 jIndex=m+n+1;
	 for k=1:Nq
		phiReac(k) = phiReac(k) + Bnm(iIndex,jIndex) * r(k)^n * Pnm(abs(m)+1,k) * exp(i*m*phi(k));
	 end
  end
end