function [dielBase,sternBase,rootDir] = readsternsrf(srfFile)

fid = fopen(srfFile);
if fid < 0
  fprintf('Could not open file %s.  Dying.\n', srfFile);
  meshBase = 0;
  return;
end
junk = fgetl(fid);
junk = fgetl(fid);
sternBase = fgetl(fid);
dielBase = fgetl(fid);
fclose(fid);

[newWord,srfFileNew] = strtok(srfFile,'/');
if length(srfFileNew) == 0
  rootDir = './';
  return;  % no directory references in srfFile.  meshBase is all
           % we need
end

if length(newWord)+length(srfFileNew)+1==length(srfFile)
  srfRootDir = sprintf('/'); % absolute filename
else 
  srfRootDir = sprintf('./'); 
end
while length(srfFileNew) > 0
  srfRootDir = sprintf('%s/%s',srfRootDir,newWord);
  [newWord,srfFileNew] = strtok(srfFileNew,'/');
end
dielBase = sprintf('%s/%s',srfRootDir,dielBase);
sternBase = sprintf('%s/%s',srfRootDir,sternBase);
rootDir = srfRootDir;