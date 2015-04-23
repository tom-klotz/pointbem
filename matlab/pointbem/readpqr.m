function pqrData = readpqr(filename,Nmax)

if nargin < 2
  Nmax = 10000;
end
xyz = zeros(Nmax,3); 
q   = zeros(Nmax,1);

fid = fopen(filename,'r');
n = 0;
while 1
  line = fgetl(fid);
  if ~ischar(line), break, end

  if length(line) > 6 ...
		  && ( strcmp(line(1:4),'ATOM') || strcmp(line(1:6),'HETATM') ...
				 )
	 
	 n = n+1;
	 [fieldName,atomNum,atomName,resName,chainID,resNum,pt,chg,radius] ...
		  = parsePQRentry(line);
	 xyz(n,:) = pt;
	 q(n) = chg;
  end
end

if n == 0
  error('Did not read any entries from PQR file!\n');
end
xyz = xyz(1:n,:);
q = q(1:n);

pqrData = struct('q', q, 'xyz', xyz);
