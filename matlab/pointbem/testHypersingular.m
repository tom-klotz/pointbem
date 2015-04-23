global kappa;
kappa = 0.01;
radius   = 2.5;
origin   = [0 0 0];
epsIn = 1;
epsOut = 80;
conv_factor = 332.112;
pqr = struct('xyz',origin,'q',1,'R',1);

minDensity = 1;
maxDensity = 20;
numDensities = 15;
densityVec = floor(logspace(log10(minDensity),log10(maxDensity),numDensities));

numHarmonics = 64;
densityVec = floor(linspace(0,2,15));

clear numpoints;
for i=1:length(densityVec)
  density = densityVec(i);
  numpoints(i) = floor(density * pi * radius^2);
  sph = makeSphereSurface(origin, radius, numpoints);
  lap = makeSurfaceToSurfaceLaplaceOperators(sph);
  [lapV,lapK,lapW] = computeLaplaceEigenvalues(radius,numHarmonics);
  lapComparison(i) = makeComparison(lap, lapV, lapK, lapW);
  yuk = makeSurfaceToSurfaceYukawaOperators(sph);
  [yukV,yukK,yukW] = computeYukawaEigenvalues(radius, kappa, numHarmonics);
  yukComparison(i) = makeComparison(yuk, yukV, yukK, yukW);
end

