function [lambdaV, lambdaK] = computeLaplaceEigenvaluesConc(sourceRadius, ...
						  destRadius, n)

highestOrderMultipole = ceil(sqrt(n));
lambdaV = zeros(highestOrderMultipole,1);
lambdaK = lambdaV;


for i=1:highestOrderMultipole
  if destRadius > sourceRadius
    ratio = sourceRadius/destRadius;
    lambdaV(1+(i-1)^2:i^2) = ratio^i * (sourceRadius/(2*(i-1)+1));
    lambdaK(1+(i-1)^2:i^2) = -2*(i-1)*ratio^i * (-1./(2.*(2*(i-1)+1.0)));
  else
    ratio = destRadius/sourceRadius;
    lambdaV(1+(i-1)^2:i^2) = ratio^(i-1) * (sourceRadius/(2*(i-1)+1));
    lambdaK(1+(i-1)^2:i^2) = 2*i*ratio^(i-1) * (-1./(2.*(2*(i-1)+1.0)));
  end
end

if destRadius > sourceRadius
  lambdaK(1) = 0;
else
  lambdaK(1) = -1;
end

lambdaV = lambdaV(1:n);
lambdaK = lambdaK(1:n);

