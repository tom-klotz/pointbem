#include <petsc.h>
#include <petsc/private/dmpleximpl.h>
//#include <slepceps.h>
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
  PetscReal startingFlops, currFlops, initializeFlops, solutionFlops;
  PetscReal estError;
  
  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);
  ierr = PetscLogDefaultBegin();CHKERRQ(ierr);
  ierr = ProcessOptions(PETSC_COMM_WORLD, &ctx);CHKERRQ(ierr);

  //get starting flops
  ierr = PetscGetFlops(&startingFlops);CHKERRQ(ierr);

  //create PQR charge distribution from PDB file or generate them if using a sphere
  if (ctx.isSphere) {
    ierr = makeSphereChargeDistribution(ctx.R, ctx.numCharges, ctx.h, PETSC_DETERMINE, &pqr);CHKERRQ(ierr);
    ierr = PQRViewFromOptions(&pqr);CHKERRQ(ierr);
  }
  else {
    PetscBool flg;
    ierr = PetscOptionsHasName(NULL, NULL, "-pqr_filename", &flg);CHKERRQ(ierr);
    if(flg) {
      ierr = PQRCreateFromPQR(PETSC_COMM_WORLD, ctx.pqrFile, &pqr);CHKERRQ(ierr);
    }
    else {
      ierr = PQRCreateFromPDB(PETSC_COMM_WORLD, ctx.pdbFile, ctx.crgFile, &pqr);CHKERRQ(ierr);
    }
  }
  
  //initialize react which stores the reaction potential at each point charge location
  ierr = VecDuplicate(pqr.q, &react);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) react, "Reaction Potential");CHKERRQ(ierr);

  //loads the mesh file
  ierr = loadSrfIntoSurfacePoints(PETSC_COMM_WORLD, ctx.srfFile, &vertNormals, &vertWeights, &panelAreas, &totalArea, &dm);CHKERRQ(ierr);

  //get flops for loading problem info
  ierr = PetscGetFlops(&currFlops);CHKERRQ(ierr);
  initializeFlops = currFlops - startingFlops;
  
  HContext params = {.alpha = ctx.alpha, .beta=ctx.beta, .gamma=ctx.gamma};

  ierr = PetscPrintf(PETSC_COMM_WORLD, "Mesh surface area: %6.6f\n", totalArea);CHKERRQ(ierr);
  PetscInt cStart, cEnd, vStart, vEnd;
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Triangular mesh has %D vertices and %D cells\n", vEnd-vStart, cEnd-cStart);CHKERRQ(ierr);
  if(ctx.usePanels == PETSC_TRUE) {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Using panel discretization...\n");CHKERRQ(ierr);
    ierr = CalculateBEMSolvationEnergy(dm, &ctx, "whocares", BEM_PANEL_MF, params, ctx.epsIn, ctx.epsOut, &pqr, panelAreas, vertNormals, react, &energy, &estError);CHKERRQ(ierr);
  }
  else {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Using point discretization...\n");CHKERRQ(ierr);
    ierr = CalculateBEMSolvationEnergy(dm, &ctx, "whocares", BEM_POINT_MF, params, ctx.epsIn, ctx.epsOut, &pqr, vertWeights, vertNormals, react, &energy, &estError);CHKERRQ(ierr);
  }

  //get flops for calculating solution including creating matrix
  ierr = PetscGetFlops(&currFlops);CHKERRQ(ierr);
  solutionFlops = currFlops - startingFlops - initializeFlops;
  
  char infoFileName[PETSC_MAX_PATH_LEN];
  PetscBool outputInfo = PETSC_FALSE;
  FILE* infoFile;
  ierr = PetscOptionsGetString(NULL, NULL, "-info_out", infoFileName, PETSC_MAX_PATH_LEN, &outputInfo);CHKERRQ(ierr);
  if(outputInfo) {
    ierr = PetscFOpen(PETSC_COMM_SELF, infoFileName, "w", &infoFile);CHKERRQ(ierr);
    if(!infoFile) {
      SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Problem opening info output file: %s", infoFileName);
    }
    //output the N and whether points or panels were used
    if(ctx.usePanels == PETSC_TRUE) {
      ierr = PetscFPrintf(PETSC_COMM_SELF, infoFile, "N = %D\nTYPE = PANELS\n", cEnd-cStart);CHKERRQ(ierr);
    }
    else {
      ierr = PetscFPrintf(PETSC_COMM_SELF, infoFile, "N=%D\nTYPE = POINTS\n", vEnd-vStart);CHKERRQ(ierr);
    }
    //output nonlinear parameters
    ierr = PetscFPrintf(PETSC_COMM_SELF, infoFile, "alpha, beta, gamma: %5.5f %5.5f %5.5f\n", params.alpha, params.beta, params.gamma);CHKERRQ(ierr);
    //output solvation energy
    ierr = PetscFPrintf(PETSC_COMM_SELF, infoFile, "solv. energy: %15.15f\n", energy);CHKERRQ(ierr);
    //output flops for each stage
    ierr = PetscFPrintf(PETSC_COMM_SELF, infoFile, "initialize FLOPS: %15.1f\nsolution FLOPS: %15.1f\n", initializeFlops, solutionFlops);CHKERRQ(ierr);
    //output residual from the last step of the BIE iteration
    ierr = PetscFPrintf(PETSC_COMM_SELF, infoFile, "final BIE residual: %5.5e\n", estError);CHKERRQ(ierr);

    //close file
    ierr = PetscFClose(PETSC_COMM_SELF, infoFile);CHKERRQ(ierr);

  }
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "We have calculated the Energy to be %.6f\n", energy);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "The total area is %6.6f\n", totalArea);CHKERRQ(ierr);

  
  ierr = PetscFinalize();
  
  return 0;
}
