loadConstants

origin = [0 0 0];
R      = 6.0;
epsIn  =  4;
epsOut = 80;
numCharges = 100;
conv_factor = 332.112;

density = 1.0;
h       = 1.0;
numPoints = ceil(4 * pi * R^2)

bem = makeBemEcfQualMatrices(surfdata, pqrdata, epsIn, epsOut);

L = bem.C * (bem.A\bem.B);
Lref = doAnalytical(R, epsIn, epsOut, pqrdata, 100); Lref = real(Lref);

E = conv_factor * 0.5 * pqrdata.q' * L * pqrdata.q;
Eref = conv_factor * 0.5 * pqrdata.q' * Lref * pqrdata.q;

fprintf('Eref = %f\nE = %f\nError = %f\nRel. error = %f\n',...
	Eref, E, Eref-E, (Eref-E)/Eref);