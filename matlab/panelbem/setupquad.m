function [ax,aw] = setupquad(N)

% sets up a 1D Gauss quadrature rule
% from Trefethen and Bau

T = zeros(N,N);
nvec = 1:N-1;
beta = 0.5 * (1. - (2.*nvec).^(-2)).^(-.5);
T = diag(beta,1)+diag(beta,-1);
[V,D]=eig(T);
x = diag(D);
w = 2*V(1,:).^2;
ax = x/2 + 0.5;
aw = w/2;

