function [jn, hn, jnprime, hnprime] = besselCai(orders, z)
% implemented from L.-W. Cai, Comp. Phys. Comm. v. 182:663-668
% (2011)

r = abs(z);
theta = atan2(imag(z),real(z));
M = floor((1.83 + 4.1*sin(theta)^0.36)*r^(0.91-0.43*sin(theta)^0.33) + 9 ...
    * (1 - sqrt(sin(theta)))); % Eq. 11 minimum M needed to get
                               % high accuracy
if M+1 < max(orders)+1
    M = max(orders)+2;
end

Mmax = floor(235 + 50 * sqrt(abs(z)));  % check against Mmax from
                                        % Eq. 12
if Mmax > 1000 
    Mmax = 1000;
end

if Mmax < M
    M = Mmax;
    if Mmax < max(orders+1)
        fprintf('error! cannot go as high as needed!\n');
    end
end

jn = besselBackwardRecurrence(0, 1e-305, M,z);
jn = jn * (sin(z)/z)/jn(1);
jnp0 = cos(z)/z - sin(z)/z^2;
jnprime = besselDerivative(jnp0,M-1,z,jn);
jn = jn(orders+1);
jnprime = jnprime(orders+1);

if imag(z) > 0
    h0 = -i*exp(i*z)/z;  % Eq 6
    h1 = (1/z -i)*h0;
    hn = besselForwardRecurrence(h0,h1,M,z);
    hnp0 = -i* (i*exp(i*z)/z - exp(i*z)/z^2); % derivative of Eq 6 with respect to argument (z)
    hnprime = besselDerivative(hnp0,M-1,z,hn);
    hn = hn(orders+1);
    hnprime = hnprime(orders+1);
else
    fprintf(['imag(z)<0 not implemented yet!\n']);
    hn = 0;
end
