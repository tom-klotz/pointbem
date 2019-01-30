#include <petsc.h>
#include <petsc/private/dmpleximpl.h>
#include "constants.h"
#include "surface.h"
#include "molecule.h"
#include "sphere.h"
#include "BEM.h"

#undef __FUNCT__
#define __FUNCT__ "ProcessOptionsSolvation"
PetscErrorCode ProcessOptionsSolvation(MPI_Comm comm, SolvationContext *ctx)
{

  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ctx->epsIn = 4; //permittivity inside molecule
  ctx->epsOut = 80; //solvent permittivity
  ctx->origin[0] = 0.0; //location of origin
  ctx->origin[1] = 0.0;
  ctx->origin[2] = 0.0;
  ctx -> h = 1.0; //charge spacing
  ctx -> density = 1.0;

  ierr = PetscOptionsBegin(comm, "", "Solvation Options", "");CHKERRQ(ierr);
  
  ierr = PetscOptionsReal("-epsilon_solute", "The dielectric coefficient of the solute", "calcSolvationEnergy", ctx->epsIn, &ctx->epsIn, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-epsilon_solvent", "The dielectric coefficient of the solvent", "calcSolvationEnergy", ctx->epsOut, &ctx->epsOut, NULL);CHKERRQ(ierr);
  //ierr = PetscOptionsString("-pdb_filename"
  
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
  DM dm;
  PQRData pqr;
  //PetscSurface msp;
  Vec panelAreas, vertWeights, vertNormals, react;
  PetscReal totalArea;
  SolvationContext ctx;
  PetscScalar energy;
  PetscErrorCode ierr;
  
  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);
  ierr = PetscLogDefaultBegin();CHKERRQ(ierr);
  ierr = ProcessOptions(PETSC_COMM_WORLD, &ctx);CHKERRQ(ierr);
  //ierr = ProcessSolvationOptions(PETSC_COMM_WORLD, &ctx);CHKERRQ(ierr);

  //create PQR charge distribution from PDB file or generate them if using a sphere
  if (ctx.isSphere) {
    ierr = makeSphereChargeDistribution(ctx.R, ctx.numCharges, ctx.h, PETSC_DETERMINE, &pqr);CHKERRQ(ierr);
    ierr = PQRViewFromOptions(&pqr);CHKERRQ(ierr);
  }
  else {
    ierr = PQRCreateFromPDB(PETSC_COMM_WORLD, ctx.pdbFile, ctx.crgFile, &pqr);CHKERRQ(ierr);
  }
  
  //initialize react which stores the reaction potential at each point charge location
  ierr = VecDuplicate(pqr.q, &react);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) react, "Reaction Potential");CHKERRQ(ierr);

  //loads the mesh file
  ierr = loadSrfIntoSurfacePoints(PETSC_COMM_WORLD, ctx.srfFile, &vertNormals, &vertWeights, &panelAreas, &totalArea, &dm);CHKERRQ(ierr);

  HContext params = {.alpha = ctx.alpha, .beta=ctx.beta, .gamma=ctx.gamma};

  ierr = CalculateBEMSolvationEnergy(dm, "whocares", BEM_POINT_MF, params, ctx.epsIn, ctx.epsOut, &pqr, vertWeights, vertNormals, react, &energy);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "We have calculated the Energy to be %.6f\n", energy);CHKERRQ(ierr);
  ierr = PetscFinalize();
  
  return 0;
}
