function [area,centroid,Z,fss,fds,fess,feds] = calcp(panel,evalpnts,directions)
% Matlab version of calcp, returns potential at evaluation point due
% to unit monopole and unit dipole uniformly distributed on a panel.
% Follows a left-hand rule (Clockwise ordered points has normal
% pointing up).
% panel -- vectors of panel vertices in rows of x,y,z (3 or 4 rows supported).
% evalpnts -- matrix of evaluation points, rows of x,y,z coordinates
% directions -- matrix of derivative directions, rows of x,y,z coordinates
% fss = the vector of potentials due to a monopole
% fds = the vector of potentials due to a panel normal dipole distribution
% area = panel area.
% centroid = panel centroid.
% Z = panel normal.
% fess = the derivative of the monopole potential at evalpnt along direction
% fess = the derivative of the dipole potential at evalpnt along direction
% 

% First check the input.
[verts, betterbethree] = size(panel);

if betterbethree ~= 3
  'wrong panel format: should be rows of x,y,z vectors!'
  return;
end

if (verts > 4) | (verts < 2) 
  'wrong panel format: panel can only have 3 or 4 vertices!'
  return;
end

if(nargin > 1) 
  [numevals, betterbethree] = size(evalpnts);

  if betterbethree ~= 3
    'wrong evaluation point format: should be rows of x,y,z vectors!'
    return;
  end
else
  numevals = 0;
end

if (nargin > 2) 
  deriv = 1;
  [numdirections, betterbethree] = size(directions);

  if betterbethree ~= 3
    'wrong direction vector format: should be rows of x,y,z vectors!'
    return;
  end

  if (numdirections ~= 0) & (numdirections ~= numevals)
    'number of direction vectors does not match number of evaluation points!'
    return;
  end
else
  deriv = 0;
end

% The PANEL SETUP!!!!!!!!!***************************************************

% Length of each side and the panel area.
for i=1:verts
  if(i < verts)
    side(i,:) = panel(i+1,:) - panel(i,:);
  else 
    side(i,:) = panel(1,:) - panel(i,:);
  end
  edgeLength(i) = norm(side(i,:));

end

% Calculate the panel coordinate system.
X = panel(3,:) - panel(1,:);
diagLength = norm(X);
if(verts == 3) 
  Y = panel(2,:) - panel(1,:); 
else 
  Y = panel(2,:) - panel(4,:); 
end

% Z-axis is normal to two diags.
Z = cross(X, Y);
area = 0.5 * norm(Z);

% Normalize panel axises. 
coord(3,:) = Z / norm(Z);
coord(1,:) = X / norm(X);
X = coord(1,:);
Z = coord(3,:);
coord(2,:) = cross(Z, X);
Y = coord(2,:);


% Determine the centroid.
vertex1 = panel(2,:) - panel(1,:);
if(verts == 4)
  vertex3 = panel(4,:) - panel(1,:);
else 
  vertex3 = panel(3,:) - panel(1,:);
end

% Project into the panel axes.
y1 = sum(vertex1 .* Y);
y3 = sum(vertex3 .* Y);
x1 = sum(vertex1 .* X);
x3 = sum(vertex3 .* X);
yc = (y1 + y3)/3.0;
xc = (diagLength + ((x1 * y1 - x3 * y3)/(y1 - y3)))/3.0;

% Compute the centroid.
centroid = panel(1,:) + xc * X + yc * Y;

% Put the corners in the newly defined coordinate system.  
for i=1:verts
  npanel(i,:) = (coord * (panel(i,:) - centroid).').';
end

% Check that panel is in the x-y panel. 
for i=1:verts
  if(abs(npanel(i,3)) > (1.0e-8 * diagLength))
    'Coordinate transform failure!!'
     npanel
     return;
  end;
end;

% Compute the contributions terms for each edge.
for i=1:verts
  if (i==verts) 
    next=1;
  else 
    next=i+1;
  end;
  ct(i) = (npanel(next,1)-npanel(i,1))/edgeLength(i);
  st(i) = (npanel(next,2)-npanel(i,2))/edgeLength(i);
end

% Done with the PANEL SETUP!!!!!!!!!************************************


% Done with Setup, now loop through the evaluation points!!
for evalindex = 1:numevals

  % Point eval pt in in new coordinate system and get the z-comp.
  point = (coord * ((evalpnts(evalindex,:) - centroid).')).';
  if(deriv == 1)
    nrm = (coord * directions(evalindex,:).').';
  end
  zn=point(3); 
  znabs=abs(zn); 
  evalDistance = norm(point);

  % Once per vertex computation
  OK = 1;
  for i=1:verts
    xc=point(1)-npanel(i,1);
    yc=point(2)-npanel(i,2);
    zc=point(3)-npanel(i,3);  
    xmxv(i)=xc; 
    ymyv(i)=yc; 
    fe(i)=xc*xc+zc*zc;
    r(i)=sqrt(yc*yc+fe(i));
    if (r(i) < (1.005*znabs)) 
      OK = 0; 
    end;
    if(deriv == 1)
      xri(i) = xmxv(i)/r(i);
      yri(i) = ymyv(i)/r(i);
    end;
  end;

  % The potential and dipole contributions are made by summing up
  % a contribution from each edge
  fs=0; 
  fd=0; 
  if(deriv == 1)
    fsx = 0; fsy = 0;
    fdx = 0; fdy = 0; fdz = 0;
  end

  for i=1:verts
    if (i==verts) 
      next=1;
    else 
      next=i+1;
    end;
    % v is the projection of the eval-i edge on the perpend to the side-i:  
    % Exploits the fact that corner points in panel coordinates. 
    v=xmxv(i)*st(i) - ymyv(i)*ct(i);

    % arg == zero if eval on next-i edge, but then v = 0. %
    arg=(r(i)+r(next)-edgeLength(i))/(r(i)+r(next)+edgeLength(i));
    if(arg == 0)
      'in calcp'
      keyboard;
    end
    fln = -log(arg);
    if (arg>0.0) 
      fs = fs + v * fln;
    end;
    if ( deriv == 1 )
      if ( arg > 0.0 ) 
        fac = (r(i)+r(next)-edgeLength(i)) * (r(i)+r(next)+edgeLength(i));
        fac = v*(edgeLength(i)+ edgeLength(i))/fac;
        fsx = fsx + (fln*st(i) - fac*(xri(i) + xri(next)));
        fsy = fsy - (fln*ct(i) + fac*(yri(i) + yri(next)));
        fdz = fdz - (fac*( 1.0/r(i) + 1.0/r(next)));
      end
    end
   
    % OK means eval not near a vertex normal, use Hess-Smith:
    if (OK == 1)
      s1=v*r(i);
      c1=znabs*(xmxv(i)*ct(i)+ymyv(i)*st(i));
      s2=v*r(next);
      c2=znabs*(xmxv(next)*ct(i)+ymyv(next)*st(i));
    else % Near a vertex normal, use Newman 
      s1=(fe(i)*st(i))-(xmxv(i)*ymyv(i)*ct(i));
      c1=znabs*r(i)*ct(i);
      s2=(fe(next)*st(i))-(xmxv(next)*ymyv(next)*ct(i));
      c2=znabs*r(next)*ct(i);
    end;
  
    s12=(s1*c2)-(s2*c1);
    c12=(c1*c2)+(s1*s2);
    val=atan2(s12, c12);
    fd=fd+val;
    if (deriv == 1) 
      u1   = xmxv(i)*ct(i) + ymyv(i)*st(i);
      u2   = xmxv(next)*ct(i)+ymyv(next)*st(i);
      if (OK == 0) % Near a vertex normal.
        rr  = r(i)*r(i);
        fh1 = xmxv(i)*ymyv(i);
        fh2 = xmxv(next)*ymyv(next);
        fac = c1/((c1*c1+s1*s1)*rr );
        if(zn < 0.0)
          fac = -1.0 * fac;
        end
        fdx = fdx + ((rr*v+fh1*u1)*fac);
        fdy = fdy - (fe(i)*u1*fac);
        rr  = r(next)*r(next);
        fac = c2/((c2*c2+s2*s2)*rr);
        if(zn < 0.0)
          fac = -1.0 * fac;
        end
        fdx = fdx - ((rr*v+fh2*u2)*fac);
        fdy = fdy + fe(next)*u2*fac;
      else 
        fac = zn/(c1*c1+s1*s1);
        fdx = fdx + (u1*v*xri(i)+r(i)*ymyv(i))*fac;
        fdy = fdy + (u1*v*yri(i)-r(i)*xmxv(i))*fac;
        fac = zn/(c2*c2+s2*s2);
        fdx = fdx - ((u2*v*xri(next)+r(next)*ymyv(next))*fac);
        fdy = fdy - ((u2*v*yri(next)-r(next)*xmxv(next))*fac);
      end
    end;
  end;
  
  if (fd<0.0) 
    fd = fd + 2*pi; 
  end;
  if (zn < 0) 
    fd=fd*(-1.0);
  end;
 
  fs=fs-zn*fd;

  if(deriv == 1 )
    fsx = fsx - zn*fdx;
    fsy = fsy - zn*fdy;
    fes = nrm(1)*fsx + nrm(2)*fsy - nrm(3)*fd;
    fed = nrm(1)*fdx + nrm(2)*fdy + nrm(3)*fdz;
  end

  % No area normalization!
  fss(evalindex) = fs;
  fds(evalindex) = fd;
  if(deriv == 1)
    fess(evalindex) = fes;
    feds(evalindex) = fed;
  end;
end;







