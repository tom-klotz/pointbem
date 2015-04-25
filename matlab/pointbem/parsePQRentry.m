function [fname,atomNum,atomName,resName,chainID,resNum,pt,chg,radius] ...
	 = parsePQRentry(line)

cpline = line;
fields = cell('');
count = 1;
while size(cpline,2) > 0
  [newTok, cpline] = strtok(cpline);
  fields(count) = cellstr(newTok);
  count = count + 1;
end
fname    = char(fields(1));
atomNum  = sscanf(char(fields(2)),'%d');
atomName = char(fields(3));
resName  = char(fields(4));

if length(fields) > 10
chainID  = char(fields(5));
resNum   = sscanf(char(fields(6)),'%d');
X        = sscanf(char(fields(7)),'%f');
Y        = sscanf(char(fields(8)),'%f');
Z        = sscanf(char(fields(9)),'%f');
chg      = sscanf(char(fields(10)),'%f');
radius   = sscanf(char(fields(11)),'%f');
else  % old style PQR files sometimes lack chainID
chainID  = '';
resNum   = sscanf(char(fields(5)),'%d');
X        = sscanf(char(fields(6)),'%f');
Y        = sscanf(char(fields(7)),'%f');
Z        = sscanf(char(fields(8)),'%f');
chg      = sscanf(char(fields(9)),'%f');
radius   = sscanf(char(fields(10)),'%f');
end  
pt = [X Y Z];