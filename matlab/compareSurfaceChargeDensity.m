addpath('./pointbem/')
addpath('./panelbem/')

loadConstants

origin = [0 0 0];
R      = 6.0;
epsIn  =  4;
epsOut = 80;
conv_factor = 332.112;

pqrdata = struct('q', 1, 'xyz', [0 0 4], 'R', R);
[gridPqrData,Ygrid,Zgrid]=  addPqrGridSpherePlane(R, struct('q',[],'xyz',[],'R',[]), 0.2);


maxOrderSphericalHarmonics = 100;
[Lref,Bnm_exact] = doAnalytical(R, epsIn, epsOut, pqrdata, maxOrderSphericalHarmonics); Lref = real(Lref);
Eref = conv_factor * 0.5 * pqrdata.q' * Lref * pqrdata.q;
analyticalReactionPotentialMap = conv_factor * computePot(Bnm_exact,gridPqrData,maxOrderSphericalHarmonics);
density = 1:2

for densityIndex = 1:length(density)
  srfFile = sprintf('pointbem/geometry/sphere_R6_vdens%d.srf',density(densityIndex)); % vdens=1:1:5
  pointSRFdata   = loadSrfIntoSurfacePoints(srfFile);
  bemPointSRF = makeBemEcfQualMatrices(pointSRFdata, pqrdata, epsIn, epsOut);
  bemPointRHS = bemPointSRF.B * pqrdata.q;
  bemPointSurfaceCharge = bemPointSRF.A \ bemPointRHS;
  bemPointReactionPotential = conv_factor * bemPointSRF.C * bemPointSurfaceCharge;
  ESRF(densityIndex) = 0.5 * pqrdata.q' * bemPointReactionPotential;



  numPoints = length(pointSRFdata.weights);
  simplesurfdata   = makeSphereSurface(origin, R, numPoints);
  bemSimple = makeBemEcfQualMatrices(simplesurfdata, pqrdata, epsIn, epsOut);
  bemSimpleRHS = bemSimple.B * pqrdata.q;
  bemSimpleSurfaceCharge = bemSimple.A \ bemSimpleRHS;
  bemSimpleReactionPotential = conv_factor * bemSimple.C * bemSimpleSurfaceCharge;
  Esimple(densityIndex) = 0.5 * pqrdata.q' * bemSimpleReactionPotential;


  panelSRFdata = loadSrfIntoPanels(srfFile);
  bemPanelSRF = makePanelBemEcfQualMatrices(panelSRFdata, pqrdata, ...
					    epsIn, epsOut);
  bemPanelRHS = bemPanelSRF.B * pqrdata.q;
  bemPanelSurfaceCharge = bemPanelSRF.A \ bemPanelRHS;
  bemPanelReactionPotential = conv_factor * bemPanelSRF.C * ...
      bemPanelSurfaceCharge;
  Epanel(densityIndex) = 0.5 * pqrdata.q' * ...
      bemPanelReactionPotential;
  
  if densityIndex==0
    bemGlobalSRF = makeBemEcfQualMatrices(pointSRFdata, gridPqrData, ...
					  epsIn, epsOut);
    bemPointReactionPotentialMap = conv_factor * bemGlobalSRF.C * bemPointSurfaceCharge;
    bemGlobalSimple = makeBemEcfQualMatrices(simplesurfdata, gridPqrData, ...
					     epsIn, epsOut);
    bemSimpleReactionPotentialMap = conv_factor * bemGlobalSimple.C * ...
	bemSimpleSurfaceCharge;
    bemGlobalPanel = makePanelBemEcfQualMatrices(panelSRFdata, gridPqrData, ...
						 epsIn, epsOut);
    bemPanelReactionPotentialMap = conv_factor * bemGlobalPanel.C * ...
	bemPanelSurfaceCharge;
  end
end
return

pots = [analyticalReactionPotentialMap; bemPointReactionPotentialMap; ...
	bemSimpleReactionPotentialMap; bemPanelReactionPotentialMap];
minPot = min(pots); maxPot = max(pots);
contourColors = minPot:(maxPot-minPot)/20:maxPot;
rPotAnalyticalGrid = reshape(analyticalReactionPotentialMap, ...
			     size(Ygrid,1), size(Ygrid,2));
rPotPointBemGrid = reshape(bemPointReactionPotentialMap,size(Ygrid, ...
						  1),size(Ygrid,2));
rPotSimpleGrid = reshape(bemSimpleReactionPotentialMap,size(Ygrid, ...
						  1),size(Ygrid,2));
rPotPanelBemGrid = reshape(bemPanelReactionPotentialMap,size(Ygrid, ...
						  1),size(Ygrid,2));


