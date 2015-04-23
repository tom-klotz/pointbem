function [Bnm,Snm,Snm2] = computeBnm_exact(b, Enm, epsIn, epsOut, Nmax)

epsHat = (epsIn-epsOut)/(0.5*(epsIn+epsOut));
for n=0:Nmax-1
  iIndex=n+1;
  for m=-n:n
	 jIndex=m+n+1;
	 Bnm(iIndex,jIndex) = ...
		  ((epsIn-epsOut)*(n+1))/(epsIn*(n*epsIn+(n+1)*epsOut)) ...
		  * (1.0/b^(2*n+1)) * Enm(iIndex,jIndex);

	 Vlambda = b/(1+2*n);
	 Klambda = -1/(2*(1+2*n));
	 Snm(iIndex,jIndex) = ...
		  Bnm(iIndex,jIndex) / Vlambda;
	 Snm2(iIndex,jIndex) = ...
		  epsHat/(1 + epsHat*Klambda) * (n+1)/b^(2*n+2)* Enm(iIndex,jIndex);
  end
end
%keyboard
