global kappa;
origin   = [0 0 0];
epsIn    = 4;
epsOut   = 80;
kappaArg = 0.0001;
r1       = 1;
r2       = 3;
pqr = struct('xyz',[0 0 0],'q',1,'R',4);

%np = [20 60 120 160 200 250 300 360 420 500 600 ];
%for i=1:length(np)
%  dielSrf = makeSphereSurface(origin, r1, np(i));%floor(density*pi*r1^2));
%  yoon = makeBemYoonDielMatrices(dielSrf, pqr, epsIn, epsOut);
%  dG(i) = 332.112 * yoon.C * (yoon.A\yoon.B) * 4*pi/2;
%end
%return

sternSrf = makeSphereSurface(origin, r2, 200);%floor(density*pi*r2^2));



stern = makeBemSternMatrices(dielSrf, sternSrf, pqr, epsIn, epsOut, kappaArg);