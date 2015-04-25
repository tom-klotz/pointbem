function bem = makeBemEcfQualMatrices(surf, pqr, epsIn, epsOut)

epsHat = (epsIn+epsOut)/(epsIn-epsOut);
surfsurfop = makeSurfaceToSurfaceLaplaceOperators(surf);
chargesurfop = makeSurfaceToChargeOperators(surf, pqr);

% the two factors below 
B = -1*    diag(surf.weights)*chargesurfop.dphidnCoul / epsIn;
C = 4*pi*  chargesurfop.slpToCharges;

A = surfsurfop.K'*diag(surf.weights) + diag(surf.weights)*epsHat/2;

bem = struct('B', B, 'A', A, 'C', C,'surfsurfop',surfsurfop,'chargesurfop',chargesurfop);