loadConstants

origin = [0 0 0];
R      = 1.0;
epsIn  =  4;
epsOut = 80;
numCharges = 10;
conv_factor = 332.112;

h       = 0.5;

pqrdata         = struct('q',1,'xyz',[0 0 0],'R',1);
Lref = doAnalytical(R, epsIn, epsOut, pqrdata, 50); Lref = real(Lref);
Eref = conv_factor * 0.5 * pqrdata.q' * Lref * pqrdata.q;
density = 1:20;
for densityIndex = 1:length(density)
  srfFile = sprintf('./geometry/sphere_R1_vdens%d.srf',density(densityIndex)); % vdens=1:1:5
  if 1
    actualSRFdata   = loadSrfIntoSurfacePoints(srfFile);
    numPointSRFUnknowns(densityIndex) = length(actualSRFdata.weights);
    bemSRF = makeBemEcfQualMatrices(actualSRFdata, pqrdata, epsIn, epsOut);
    LSRF = bemSRF.C * (bemSRF.A\bemSRF.B);
    ESRF(densityIndex) = conv_factor * 0.5 * pqrdata.q' * LSRF * pqrdata.q;
  end
  
  if 1 
    scaleFactor(densityIndex) = (4 * pi * R^2)/sum(actualSRFdata.weights);
    newSRFdata = actualSRFdata;
    newSRFdata.weights = scaleFactor(densityIndex) * actualSRFdata.weights;
    bemScaledSRF = makeBemEcfQualMatrices(newSRFdata, pqrdata, epsIn, epsOut);
    LScaledSRF = bemScaledSRF.C * (bemScaledSRF.A\bemScaledSRF.B);
    EScaledSRF(densityIndex) = conv_factor * 0.5 * pqrdata.q' * LScaledSRF * pqrdata.q;
  end

  if 1
    numSimpleUnknowns(densityIndex) = ceil(4 * pi * density(densityIndex)*R^2);
    simplesurfdata   = makeSphereSurface(origin, R, numSimpleUnknowns(densityIndex));
    bemSimple = makeBemEcfQualMatrices(simplesurfdata, pqrdata, epsIn, epsOut);
    LSimple = bemSimple.C * (bemSimple.A\bemSimple.B);
    ESimple(densityIndex) = conv_factor * 0.5 * pqrdata.q' * LSimple * pqrdata.q;
  end
end

figure;
set(gca,'fontsize',16);
loglog(numPointSRFUnknowns,abs(Eref-ESRF),'b','linewidth',2);
hold on
loglog(numPointSRFUnknowns,abs(Eref-EScaledSRF),'r','linewidth',2);
loglog(numSimpleUnknowns,abs(Eref-ESimple),'k','linewidth',2);


%fprintf('Eref = %f\nESRF = %f\nError = %f\nRel. error = %f\n',...
%	Eref, ESRF, Eref-ESRF, (Eref-ESRF)/Eref);
%fprintf('Eref = %f\nESimple = %f\nError = %f\nRel. error = %f\n',...
%	Eref, ESimple, Eref-ESimple, (Eref-ESimple)/Eref);

