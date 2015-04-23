function [R, t, Bg] = makePanelTransform(panel)


% This does everything with points as column vectors. just FYI

op1 = panel(1,:)';
op2 = panel(2,:)';
op3 = panel(3,:)';
t = op1;

p1 = op1-t;
p2 = op2-t;
p3 = op3-t;

x = p2-p1; x=x/norm(x);
y = p3-p1; y=y-x*(x'*y);y=y/norm(y);
z = cross(x,y);

R = [x'; y'; z']';
t=t';

% at this point, if you called this function with panel p,
% and set ptrans=(p-repmat(t,[3,1]))*R, then
% ptrans  = [ 0 0 0;
%            x2 0 0;
%    	     x3 y3 0;];
% and you recover p with
% p3 = ptrans*R'+repmat(t,[3,1])

% we also need the affine transformation to the unit triangle.  this
% is just ptrans(2:3,1:2)!

paneltransform=(panel-repmat(t,[3,1]))*R;
Bg = paneltransform(2:3,1:2);

% to compute global coord from ref triangle pt [eta ksi]:
% globalCoord = [[eta ksi]*Bg 0]*R' +t

