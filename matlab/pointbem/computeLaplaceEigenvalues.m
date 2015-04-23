function [lambdaV, lambdaK, lambdaW] = computeLaplaceEigenvalues(radius, n)

highestOrderMultipole = ceil(sqrt(n));
lambdaV = zeros(highestOrderMultipole,1);
lambdaK = lambdaV;
lambdaW = lambdaV;
for i=1:highestOrderMultipole
  lambdaV(1+(i-1)^2:i^2) = radius/(2*(i-1)+1);
  lambdaK(1+(i-1)^2:i^2) = -1./(2.*(2*(i-1)+1.0));
  lambdaW(1+(i-1)^2:i^2) = ((i-1)*(i))/(radius * (2*(i-1)+1));
end


lambdaV = lambdaV(1:n);
lambdaK = lambdaK(1:n);
lambdaW = lambdaW(1:n);
