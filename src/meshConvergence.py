#!/usr/bin/env python
import utils

def plotConv(args):
  from pylab import legend, semilogx, loglog, plot, savefig, show, title, xlabel, ylabel

  output = utils.parseOutput(args)
  Nc_SRF = output[:,2]
  Np_SRF = output[:,0]
  E_SRF  = output[:,4]
  E_PAN  = output[:,6]
  if args.type == 'sphere':
    E_ref = output[:,3]
  else:
    E_ref = utils.rich(E_PAN[-2], E_PAN[-1], 2, 1)

  Err_SRF    = E_ref - E_SRF
  Err_PAN    = E_ref - E_PAN
  RelErr_SRF = Err_SRF / E_ref
  RelErr_PAN = Err_PAN / E_ref

  RelErr_Sqrt = 7e-1*(1.0/Nc_SRF)**0.5
  RelErr_Lin  = 5e0*(1.0/Nc_SRF)

  loglog(Np_SRF, abs(RelErr_SRF), Nc_SRF, abs(RelErr_PAN), Np_SRF, RelErr_Sqrt, Nc_SRF, RelErr_Lin)
  if args.type == 'sphere':
    title(args.type.upper())
  else:
    title('Residue: '+args.type.upper())
  xlabel('Number of Dof')
  ylabel('Relative Error')
  legend(['SRF', 'PAN', 'Sqrt', 'Linear'], 'upper right', shadow = True)
  if not args.output is None:
    savefig(args.output+'_RelErr.png')
  else:
    show()

  semilogx(Np_SRF, E_SRF, Nc_SRF, E_PAN, Nc_SRF, [E_ref for i in range(len(Nc_SRF))])
  if args.type == 'sphere':
    title(args.type.upper())
  else:
    title('Residue: '+args.type.upper())
  xlabel('Number of Dof')
  ylabel('Solvation Energy (kcal/mol)')
  legend(['SRF', 'PAN', 'REF'], 'lower right', shadow = True)
  if not args.output is None:
    savefig(args.output+'_SolvEng.png')
  else:
    show()

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(description     = 'Work Precision Diagrams for BEM',
                                   epilog          = 'For more information, visit http://www.bitbucket.org/knepely/pointbem-petsc',
                                   formatter_class = argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--type', default='sphere', help='The run type, e.g. sphere, arg, asp, ...')
  parser.add_argument('--output', help='Basename for output plots')

  args = parser.parse_args()
  print(args)
  plotConv(args)
