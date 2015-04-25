global kappa;
kappa = 0.01;
radius   = 2.5;
origin   = [0 0 0];
epsIn = 1;
epsOut = 80;
conv_factor = 332.112;
pqr = struct('xyz',origin,'q',1,'R',1);
dG_exact = 0.5 * conv_factor * (1/epsOut - 1/epsIn) * 1/radius;

minDensity = 1;
maxDensity = 20;
numDensities = 15;
densityVec = floor(logspace(log10(minDensity),log10(maxDensity),numDensities));

for i=1:numDensities
  density = densityVec(i);
  numpoints(i) = floor(density * pi * radius^2);
  sph = makeSphereSurface(origin, radius, numpoints(i));
  yoon = makeBemYoonDielMatrices(sph, pqr, epsIn, epsOut);
  dG_solv(i) = 0.5 * conv_factor * yoon.C * (yoon.A\yoon.B);
end

figure;
set(gca,'fontsize',16);
loglog(numpoints, abs(dG_solv-dG_exact),'s--','linewidth',2);
hold on;
loglog([numpoints(1) 10*numpoints(1)], ...
       abs(dG_solv(1)-dG_exact)*2* [1 1/sqrt(10)],'k-.','linewidth',2);
xlabel('N_{points}');
ylabel('Error (kcal/mol)');
