function pqr = makeGridCharge(L, h)

x = -L/2:h:L/2;
[X,Y,Z] = meshgrid(x, x, x);
xyz = [reshape(X,numel(X),1) reshape(Y,numel(Y),1) reshape(Z, numel(Z),1)];
pqr = struct('q',zeros(numel(X),1),'xyz',xyz,'r',zeros(numel(X),1));
