function surfdata = loadSrfIntoPanels(srf)

[dielBase,sternBase,rootDir] = readsternsrf(srf);
surfdata = loadAndProcessMesh(dielBase, 1, 1);
