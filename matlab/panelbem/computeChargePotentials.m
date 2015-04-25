function [slp,dlp] = computeChargePotentials(meshData,chargeLocations)

np = size(meshData.face,1);
ncharges = size(chargeLocations,1);

slp = zeros(ncharges,np);
dlp = slp;

for i=1:np
  panel = [meshData.X(:,i) meshData.Y(:,i) meshData.Z(:,i)];
  [area,collocpt,Z,fss,fds] = calcp(panel,chargeLocations);
  slp(:,i) = fss'/4/pi;
  dlp(:,i) = fds'/4/pi;
end
