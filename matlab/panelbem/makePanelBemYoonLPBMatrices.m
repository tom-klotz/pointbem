function bem = makePanelBemYoonLPBMatrices(surf, pqr, epsIn, ...
					    epsOut, kappaArg)
global kappa;
kappa = kappaArg;

npDiel = length(surf.areas);
nCharges = length(pqr.q);

dielChargeOp = makeSurfaceToChargePanelOperators(surf,pqr);

Idiel = eye(npDiel);

dielDielOp = makeSurfaceToSurfaceLaplacePanelOperators(surf);
dielDielYukawa  = makeSurfaceToSurfaceYukawaPanelOperators(surf,kappaArg);

A11 = (Idiel/2 + dielDielOp.K);
A12 = -dielDielOp.V;
A21 = (Idiel/2 - dielDielYukawa.K);
A22_base = dielDielYukawa.V;
A = [A11 A12; 
     A21 (epsIn/epsOut)*A22_base;];

B = [dielChargeOp.phiCoul/epsIn;
     zeros(npDiel,nCharges)];

C = 4 * pi * [-dielChargeOp.dlpToCharges dielChargeOp.slpToCharges];

bem = struct('B', B, 'A', A, 'C', C,...
	     'dielChargeOp',dielChargeOp,...
	     'dielDielOp',dielDielOp,...
	     'dielDielYukawa',dielDielYukawa,...
	     'A11',A11,'A12',A12,...
	     'A21',A21,'A22_base',A22_base);