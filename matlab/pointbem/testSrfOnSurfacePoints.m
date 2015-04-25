loadConstants

origin = [0 0 0];
R      = 6.0;
srfFile = './geometry/sphere_R6_vdens1.srf'; % vdens=1:1:5
epsIn  =  4;
epsOut = 80;
numCharges = 100;
conv_factor = 332.112;

density = 1.0;
h       = 1.0;
numPoints = ceil(4 * pi * R^2)

pqrdata         = makeSphereChargeDistribution(R, numCharges, h); 

actualSRFdata   = loadSrfIntoSurfacePoints(srfFile);
bemSRF = makeBemEcfQualMatrices(actualSRFdata, pqrdata, epsIn, epsOut);

simplesurfdata   = makeSphereSurface(origin, R, numPoints);
bemSimple = makeBemEcfQualMatrices(simplesurfdata, pqrdata, epsIn, epsOut);

LSRF = bemSRF.C * (bemSRF.A\bemSRF.B);
LSimple = bemSimple.C * (bemSimple.A\bemSimple.B);

Lref = doAnalytical(R, epsIn, epsOut, pqrdata, 100); Lref = real(Lref);
Eref = conv_factor * 0.5 * pqrdata.q' * Lref * pqrdata.q;

ESRF = conv_factor * 0.5 * pqrdata.q' * LSRF * pqrdata.q;
ESimple = conv_factor * 0.5 * pqrdata.q' * LSimple * pqrdata.q;

fprintf('Eref = %f\nESRF = %f\nError = %f\nRel. error = %f\n',...
	Eref, ESRF, Eref-ESRF, (Eref-ESRF)/Eref);
fprintf('Eref = %f\nESimple = %f\nError = %f\nRel. error = %f\n',...
	Eref, ESimple, Eref-ESimple, (Eref-ESimple)/Eref);

