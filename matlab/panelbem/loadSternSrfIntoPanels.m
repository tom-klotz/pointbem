function surfdata = loadSternSrfIntoPanels(srf)


[dielBase,sternBase,rootDir] = readsternsrf(srf);

dielBndy(1) = loadAndProcessMesh(dielBase, 1, 1); 
sternBndy(1) = loadAndProcessMesh(sternBase, 1, 1);

surfdata = struct('dielBndy',dielBndy,'sternBndy',sternBndy,...
		  'numDielPanels',[length(dielBndy(1).areas)],...
		  'numSternPanels',[length(sternBndy(1).areas)]);

