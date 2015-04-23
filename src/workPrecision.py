#!/usr/bin/env python
import utils

def plot(args):
  from pylab import legend, semilogx, loglog, savefig, show, title, xlabel, ylabel

  output = utils.parseOutput(args)
  Nc_SRF = output[:,2]
  Np_SRF = output[:,0]
  Np_PNT = output[:,1]
  E_SRF  = output[:,4]
  E_PNT  = output[:,5]
  E_PAN  = output[:,6]
  if args.type == 'sphere':
    E_ref = output[:,3]
  else:
    E_ref = utils.rich(E_PAN[-2], E_PAN[-1], 2, 1)
  Flops_SRF     = output[:,7]
  FlopsS2S_SRF  = output[:,8]
  Flops_PAN     = output[:,9]
  FlopsS2S_PAN  = output[:,10]
  Err_SRF       = E_ref - E_SRF
  Err_PNT       = E_ref - E_PNT
  Err_PAN       = E_ref - E_PAN

  doTotal = True
  if doTotal:
    semilogx(Flops_SRF, abs(Err_SRF), 'r--', Flops_PAN, abs(Err_PAN), 'b--', FlopsS2S_SRF, abs(Err_SRF), 'r', FlopsS2S_PAN, abs(Err_PAN), 'b')
    legend(['SRF', 'PAN', 'SRF S2S', 'PAN S2S'], 'upper right', shadow = True)
  else:
    semilogx(FlopsS2S_SRF, abs(Err_SRF), 'r', FlopsS2S_PAN, abs(Err_PAN), 'b')
    legend(['SRF S2S', 'PAN S2S'], 'upper right', shadow = True)
  if args.type == 'sphere':
    title(args.type.upper())
  else:
    title('Residue: '+args.type.upper())
  xlabel('Number of Flops')
  ylabel('Absolute Error (kcal/mol)')
  if not args.output is None:
    savefig(args.output+'.png')
  else:
    show()
  return

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(description     = 'Work Precision Diagrams for BEM',
                                   epilog          = 'For more information, visit http://www.bitbucket.org/knepely/pointbem-petsc',
                                   formatter_class = argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--type', default='sphere', help='The run type, e.g. sphere, arg, asp, ...')
  parser.add_argument('--output', help='Basename for output plots')

  args = parser.parse_args()
  print(args)
  plot(args)
