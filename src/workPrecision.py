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
  fsurfRE   = r'Flops_Surf\s+= (?P<flops>[0-9]\.[0-9]+e\+[0-9][0-9])\s+Flops_S2S_Surf\s+= (?P<flops_s2s>[0-9]\.[0-9]+e\+[0-9][0-9])'
  fpanelRE  = r'Flops_Panel\s+= (?P<flops>[0-9]\.[0-9]+e\+[0-9][0-9])\s+Flops_S2S_Panel\s+= (?P<flops_s2s>[0-9]\.[0-9]+e\+[0-9][0-9])'
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
      m = re.match(fsurfRE, line)
      if m:
        fsurf    = m.group('flops')
        fsurfS2S = m.group('flops_s2s')
        continue
      m = re.match(fpanelRE, line)
      if m:
        fpanel    = m.group('flops')
        fpanelS2S = m.group('flops_s2s')
        continue
  return int(verts), int(vertsSim), int(cells), float(eref), float(esurf), float(esimple), float(epanel), float(fsurf), float(fsurfS2S), float(fpanel), float(fpanelS2S)

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
Flops_SRF     = output[:,7]
FlopsS2S_SRF  = output[:,8]
Flops_PAN     = output[:,9]
FlopsS2S_PAN  = output[:,10]
Err_SRF       = E_ref - E_SRF
Err_PNT       = E_ref - E_PNT
Err_PAN       = E_ref - E_PAN

from pylab import legend, semilogx, loglog, savefig, show, title, xlabel, ylabel

semilogx(Flops_SRF, abs(Err_SRF), Flops_PAN, abs(Err_PAN), FlopsS2S_SRF, abs(Err_SRF), FlopsS2S_PAN, abs(Err_PAN))
title('Work-Precision Diagram')
xlabel('Number of Flops')
ylabel('Absolute Error (kcal/mol)')
legend(['SRF', 'PAN', 'SRF S2S', 'PAN S2S'], 'upper right', shadow = True)
if len(sys.argv) > 1:
  savefig(sys.argv[1])
else:
  show()
