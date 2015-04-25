function bem = makePanelBemEcfQualMatrices(surf, pqr, epsIn, epsOut)

epsHat = (epsIn+epsOut)/(epsIn-epsOut);
surfsurfop = makeSurfaceToSurfaceLaplacePanelOperators(surf);
chargesurfop = makeSurfaceToChargePanelOperators(surf, pqr);

% the two factors below 
B = -1*    chargesurfop.dlpToCharges' / epsIn;
C = 4*pi*  chargesurfop.slpToCharges;

A = surfsurfop.K'*diag(surf.areas) + diag(surf.areas)*epsHat/2;

bem = struct('B', B, 'A', A, 'C', C,'surfsurfop',surfsurfop,'chargesurfop',chargesurfop);

