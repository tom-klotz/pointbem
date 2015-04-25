function [lambdaV, lambdaK, lambdaW] = computeYukawaEigenvalues(radius, kappa, n)

highestOrderMultipole = ceil(sqrt(n));
lambdaV = zeros(highestOrderMultipole,1);
lambdaK = lambdaV;
lambdaW = lambdaV;

[jn3,hn3,jnprime3,hnprime3]=besselCai(0:highestOrderMultipole,i*kappa*radius);
vlist3 = real(i*(i*kappa*radius^2*jn3.*hn3));
klist3 = real(i*((i*kappa)^2*radius^2/2)*(jnprime3.*hn3 + hnprime3.*jn3));
wlist3 = real(-i*((i*kappa)^3*radius^2)*(jnprime3.*hnprime3));

for j=1:highestOrderMultipole
  lambdaV(1+(j-1)^2:j^2) = vlist3(j);
  lambdaK(1+(j-1)^2:j^2) = klist3(j);
  lambdaW(1+(j-1)^2:j^2) = wlist3(j);
end
lambdaV = lambdaV(1:n);
lambdaK = lambdaK(1:n);
lambdaW = lambdaW(1:n);
%fprintf('in computeYukawa: x = %f\n',kappa*radius);
