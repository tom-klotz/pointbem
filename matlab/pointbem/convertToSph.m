function sphpnt = convertToSph(pnt)

sphpnt = zeros(size(pnt,1),3);
for i=1:size(pnt,1)
  rho = norm(pnt(i,:));
  theta = atan2(pnt(i,2), pnt(i,1));
  phi = 0;
  if abs(rho) > 0
	 phi = acos(pnt(i,3)/rho);
  end
  sphpnt(i,:) = [rho theta phi];
end