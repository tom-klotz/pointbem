function bem = makePanelBemSternMatrices(surf, pqr, epsIn, epsOut, ...
					 kappaArg)
global kappa;
kappa = kappaArg;

npDiel = length(surf.dielBndy(1).areas);
npStern = length(surf.sternBndy(1).areas);
npTotal = npDiel+npStern;

nCharges = length(pqr.q);

dielChargeOp = makeSurfaceToChargePanelOperators(surf.dielBndy(1), ...
						 pqr);
Idiel = eye(npDiel);
Istern = eye(npStern);

% src = diel, target = diel 
dielDielOp = makeSurfaceToSurfaceLaplacePanelOperators(surf.dielBndy(1));

% src = diel, target = stern
dielToSternOp = makeSurfaceToSurfaceLaplacePanelOperators(surf.dielBndy(1), ...
						  surf.sternBndy(1));

% src = stern, target = diel
sternToDielOp = ...
    makeSurfaceToSurfaceLaplacePanelOperators(surf.sternBndy(1), ...
					      surf.dielBndy(1));

% LAPLACE kernel for src = stern, target = stern
sternSternLaplace = makeSurfaceToSurfaceLaplacePanelOperators(surf.sternBndy(1));

% YUKAWA kernel for src = stern, target = stern
if kappa > 1e-6
  sternSternYukawa = ...
      makeSurfaceToSurfaceYukawaPanelOperators(surf.sternBndy(1), kappa);
else
  sternSternYukawa = ...
      makeSurfaceToSurfaceLaplacePanelOperators(surf.sternBndy(1));
end

Z1 = zeros(npDiel, npStern);
Z2 = zeros(npStern, npDiel); 
% could just use Z1' but this is clearer

A11 = (Idiel/2 + dielDielOp.K);
A12 = -dielDielOp.V;
A13 = Z1;
A14 = Z1;

A21 = (Idiel/2 - dielDielOp.K);
A22_base = dielDielOp.V;
A23 = sternToDielOp.K;
A24 = -sternToDielOp.V;

A31 = -dielToSternOp.K;
A32_base = dielToSternOp.V;
A33 = (Istern/2+sternSternLaplace.K);
A34 = -sternSternLaplace.V;

A41 = Z2;
A42 = Z2; 
A43 = (Istern/2 - sternSternYukawa.K);
A44_base = sternSternYukawa.V;

A = [A11 A12 A13 A14;
     A21 (epsIn/epsOut)*A22_base A23 A24;
     A31 (epsIn/epsOut)*A32_base A33 A34;
     A41 A42 A43 (epsOut/epsOut)*A44_base];

B = [dielChargeOp.phiCoul/epsIn;
     zeros(npDiel,nCharges);
     zeros(npStern,nCharges);
     zeros(npStern,nCharges);];

C = 4*pi*[-dielChargeOp.dlpToCharges dielChargeOp.slpToCharges ...
     zeros(nCharges,npStern) zeros(nCharges, npStern)];

bem = struct('B', B, 'A', A, 'C', C,...
	     'dielChargeOp',dielChargeOp,...
	     'dielDielOp',dielDielOp,...
	     'dielToSternOp',dielToSternOp,...
	     'sternToDielOp',sternToDielOp,...
	     'sternSternLaplace',sternSternLaplace,...
	     'sternSternYukawa',sternSternYukawa,...
	     'A11',A11,'A12',A12,'A13',A13,'A14',A14,...
	     'A21',A21,'A22_base',A22_base,'A23',A23,'A24',A24,...
	     'A31',A31,'A32_base',A32_base,'A33',A33,'A34',A34,...
	     'A41',A41,'A42',A42,'A43',A43,'A44_base',A44_base);
