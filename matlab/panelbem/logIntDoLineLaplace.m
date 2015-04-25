function [phi,phiDouble] = logIntDoLineLaplace(kappa,z,x,y1,y2,foo)

global quadrule_order quadrule_x quadrule_w;

phi = 0;
phiDouble = 0;
%fprintf('logInt: %f, %f, %f, %f, %f\n',...
%		  kappa,z,x,y1,y2);
%return;

theta1 = atan2(y1,x);
theta2 = atan2(y2,x);

dtheta = theta2 - theta1;
absZ = abs(z);
signZ = sign(z);
for i=1:quadrule_order
  theta = theta1 + dtheta*quadrule_x(i);
  weight =         dtheta*quadrule_w(i);
  Rtheta = x * sec(theta);
  r = sqrt(Rtheta^2+z^2);
  functionAtTheta(i) = r - absZ;
  doubleFunctionAtTheta(i) = z/r - signZ;
  phi = phi + weight * functionAtTheta(i);
  phiDouble = phiDouble + weight*doubleFunctionAtTheta(i);
end


if nargin > 5
%  fprintf('th1=%f  th2=%f\n',theta1,theta2);
%  keyboard;
end
