def rich(coarse, fine, k, n):
  # Richardson extrapolation: http://en.wikipedia.org/wiki/Richardson_extrapolation
  # R(h, k) = \frac{k^n A(h) - A(k h)}{k^n - 1}
  return (k**n * fine - coarse)/(k**n - 1)

def parseOutputFile(filename):
  import re
  vertsSim = -1
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

def parseOutput(args):
  import numpy as np
  import glob, os, re

  fnRE   = r'output/testSrfOnSurfacePoints_[a-zA-Z]+_(?P<num>\d+)\.out'
  files  = glob.glob(os.path.join('output', 'testSrfOnSurfacePoints_%s_*.out' % args.type))
  files.sort(key = lambda name: int(re.search(fnRE, name).group('num')))
  output = [parseOutputFile(fn) for fn in files]
  return np.array(output)
