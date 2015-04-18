#!/usr/bin/env python
import numpy as np
import sys

def parseOutput(filename):
  import re
  eref = esurf = esimple = epanel = -1
  sizRE     = r'SRF (?P<verts>\d+) vertices (?P<cells>\d+) cells'
  sizSimRE  = r'Simple (?P<verts>\d+) vertices'
  esurfRE   = r'Eref = (?P<ref>[+-]?[0-9]+\.[0-9]+) ESurf\s+= (?P<energy>[+-]?[0-9]+\.[0-9]+) Error = (?P<error>[+-]?[0-9]+\.[0-9]+) Rel. error = (?P<relerror>[+-]?[0-9]+\.[0-9]+)'
  esimpleRE = r'Eref = (?P<ref>[+-]?[0-9]+\.[0-9]+) ESimple\s+= (?P<energy>[+-]?[0-9]+\.[0-9]+) Error = (?P<error>[+-]?[0-9]+\.[0-9]+) Rel. error = (?P<relerror>[+-]?[0-9]+\.[0-9]+)'
  epanelRE  = r'Eref = (?P<ref>[+-]?[0-9]+\.[0-9]+) EPanel\s+= (?P<energy>[+-]?[0-9]+\.[0-9]+) Error = (?P<error>[+-]?[0-9]+\.[0-9]+) Rel. error = (?P<relerror>[+-]?[0-9]+\.[0-9]+)'
  with file(filename) as f:
    for line in f.readlines():
      m = re.match(sizRE, line)
      if m:
        verts = m.group('verts')
        cells = m.group('cells')
        continue
      m = re.match(sizSimRE, line)
      if m:
        vertsSim = m.group('verts')
        continue
      m = re.match(esurfRE, line)
      if m:
        eref  = m.group('ref')
        esurf = m.group('energy')
        continue
      m = re.match(esimpleRE, line)
      if m:
        esimple = m.group('energy')
        continue
      m = re.match(epanelRE, line)
      if m:
        epanel = m.group('energy')
        continue
  return int(verts), int(vertsSim), int(cells), float(eref), float(esurf), float(esimple), float(epanel)

basename = 'output/testSrfOnSurfacePoints_%d.out'
output   = [parseOutput(basename % (tnum)) for tnum in range(1, 9)]
output   = np.array(output)

Nc_SRF = output[:,2]
Np_SRF = output[:,0]
Np_PNT = output[:,1]
E_SRF  = output[:,4]
E_PNT  = output[:,5]
E_PAN  = output[:,6]
E_ref  = output[:,3]

Err_SRF    = E_ref - E_SRF
Err_PNT    = E_ref - E_PNT
Err_PAN    = E_ref - E_PAN
RelErr_SRF = Err_SRF / E_ref
RelErr_PNT = Err_PNT / E_ref
RelErr_PAN = Err_PAN / E_ref

RelErr_Sqrt = 7e-1*(1.0/Nc_SRF)**0.5
RelErr_Lin  = 5e0*(1.0/Nc_SRF)
RelErr_Quad = 3e4*(1.0/Nc_SRF)**2

from pylab import legend, loglog, savefig, show, title, xlabel, ylabel

loglog(Np_SRF, abs(RelErr_SRF), Np_PNT, abs(RelErr_PNT), Nc_SRF, abs(RelErr_PAN), Np_SRF, RelErr_Sqrt, Nc_SRF, RelErr_Lin)
title('Mesh Convergence')
xlabel('Number of Dof')
ylabel('Relative Error')
legend(['SRF', 'PNT', 'PAN', 'Sqrt', 'Linear'], 'upper left', shadow = True)
if len(sys.argv) > 1:
  savefig(sys.argv[1])
else:
  show()
