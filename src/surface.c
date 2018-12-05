#include <petsc.h>
#include "constants.h"
#include "surface.h"

#undef __FUNCT__
#define __FUNCT__ "DMPlexCreateBardhanFromFile"
/*@C
  DMPlexCreateBardhanFromFile - Create a DMPlex mesh from a Bardhan surface mesh file

+ comm        - The MPI communicator
. filename    - Name of the Fluent mesh file
- interpolate - Create faces and edges in the mesh

  Output Parameters:
+ n  - The vertex normals
- dm - The DM object representing the mesh

  Level: beginner

.seealso: DMPlexCreateFromFile(), DMPlexCreateFluent(), DMPlexCreateExodus(), DMPlexCreate()
@*/
PetscErrorCode DMPlexCreateBardhanFromFile(MPI_Comm comm, const char filename[], PetscBool interpolate, Vec *n, DM *dm)
{
  PetscViewer     viewerVert, viewerFace;
  char            vertname[PETSC_MAX_PATH_LEN];
  char            facename[PETSC_MAX_PATH_LEN];
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = PetscViewerCreate(comm, &viewerVert);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewerVert, PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewerVert, FILE_MODE_READ);CHKERRQ(ierr);
  ierr = PetscStrcpy(vertname, filename);CHKERRQ(ierr);
  ierr = PetscStrcat(vertname, ".vert");CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewerVert, vertname);CHKERRQ(ierr);
  ierr = PetscViewerCreate(comm, &viewerFace);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewerFace, PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewerFace, FILE_MODE_READ);CHKERRQ(ierr);
  ierr = PetscStrcpy(facename, filename);CHKERRQ(ierr);
  ierr = PetscStrcat(facename, ".face");CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewerFace, facename);CHKERRQ(ierr);
  ierr = DMPlexCreateBardhan(comm, viewerVert, viewerFace, interpolate, n, dm);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewerVert);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewerFace);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPlexCreateBardhan"
/*@C
  DMPlexCreateBardhan - Create a DMPlex mesh from a Bardhan surface mesh file.

  Collective on comm

  Input Parameters:
+ comm  - The MPI communicator
. viewerVert - The Viewer associated with a vertex file
. viewerFace - The Viewer associated with a face file
- interpolate - Create faces and edges in the mesh

  Output Parameters:
+ n  - The vertex normals
- dm - The DM object representing the mesh

  Level: beginner

.keywords: surface mesh
.seealso: DMPLEX, DMCreate()
@*/
PetscErrorCode DMPlexCreateBardhan(MPI_Comm comm, PetscViewer viewerVert, PetscViewer viewerFace, PetscBool interpolate, Vec *n, DM *dm)
{
  PetscSection   coordSection;
  Vec            coordinates;
  PetscScalar   *cds;
  PetscReal     *coords, *normals;
  PetscInt      *cells;
  PetscInt       dim = 2, dimEmbed = 3, numCells = 0, numVertices = 0, numCellVertices = 3, coordSize, c, v, d;
  PetscMPIInt    rank;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  if (!rank) {
    PetscInt maxSize = 1000, cnt, cnt2, tot;

    ierr = PetscMalloc2(maxSize, &coords, maxSize, &normals);CHKERRQ(ierr);
    cnt  = 3;
    tot  = 0;
    while (cnt == 3) {
      PetscInt dummy[3];

      if (tot + cnt > maxSize) {
        PetscReal *tmpcoords, *tmpnormals;

        ierr = PetscMalloc2(maxSize*2, &tmpcoords, maxSize*2, &tmpnormals);CHKERRQ(ierr);
        ierr = PetscMemcpy(tmpcoords,  coords,  maxSize * sizeof(PetscReal));CHKERRQ(ierr);
        ierr = PetscMemcpy(tmpnormals, normals, maxSize * sizeof(PetscReal));CHKERRQ(ierr);
        ierr = PetscFree2(coords, normals);CHKERRQ(ierr);
        coords  = tmpcoords;
        normals = tmpnormals;
        maxSize *= 2;
      }
      ierr = PetscViewerRead(viewerVert, &coords[tot],  3, &cnt, PETSC_DOUBLE);CHKERRQ(ierr);
      ierr = PetscViewerRead(viewerVert, &normals[tot], 3, &cnt, PETSC_DOUBLE);CHKERRQ(ierr);
      ierr = PetscViewerRead(viewerVert, dummy,         3, &cnt, PETSC_INT);CHKERRQ(ierr);
      tot += cnt;
      if (cnt == 3) ++numVertices;
    }

    cnt  = 3;
    tot  = 0;
    ierr = PetscMalloc1(maxSize, &cells);CHKERRQ(ierr);
    while (cnt == 3) {
      PetscInt dummy[3];

      if (tot + cnt > maxSize) {
        PetscInt *tmpcells;

        ierr = PetscMalloc1(maxSize*2, &tmpcells);CHKERRQ(ierr);
        ierr = PetscMemcpy(tmpcells, cells, maxSize * sizeof(PetscInt));CHKERRQ(ierr);
        ierr = PetscFree(cells);CHKERRQ(ierr);
        cells = tmpcells;
        maxSize *= 2;
      }
      ierr = PetscViewerRead(viewerFace, &cells[tot], 3, &cnt,  PETSC_INT);CHKERRQ(ierr);
      ierr = PetscViewerRead(viewerFace, dummy,       2, &cnt2, PETSC_INT);CHKERRQ(ierr);
      tot += cnt;
      if (cnt == 3) ++numCells;
    }
    for (c = 0; c < tot; ++c) cells[c] += numCells-1;
  }

  /* Allocate cell-vertex mesh */
  ierr = DMCreate(comm, dm);CHKERRQ(ierr);
  ierr = DMSetType(*dm, DMPLEX);CHKERRQ(ierr);
  ierr = DMSetDimension(*dm, dim);CHKERRQ(ierr);
  ierr = DMPlexSetChart(*dm, 0, numCells + numVertices);CHKERRQ(ierr);
  if (numCells < 0) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid number of cells %d in Bardhan file", numCells);
  for (c = 0; c < numCells; ++c) {ierr = DMPlexSetConeSize(*dm, c, numCellVertices);CHKERRQ(ierr);}
  ierr = DMSetUp(*dm);CHKERRQ(ierr);
  for (c = 0; c < numCells; ++c) {
    ierr = DMPlexSetCone(*dm, c, &cells[c*numCellVertices]);CHKERRQ(ierr);
  }
  ierr = DMPlexSymmetrize(*dm);CHKERRQ(ierr);
  ierr = DMPlexStratify(*dm);CHKERRQ(ierr);
  if (interpolate) {
    DM idm = NULL;

    ierr = DMPlexInterpolate(*dm, &idm);CHKERRQ(ierr);
    ierr = DMDestroy(dm);CHKERRQ(ierr);
    *dm  = idm;
  }

  /* TODO: What do Jay's vertex checks do? */

  /* Read coordinates */
  ierr = DMSetCoordinateDim(*dm, dimEmbed);CHKERRQ(ierr);
  ierr = DMGetCoordinateSection(*dm, &coordSection);CHKERRQ(ierr);
  ierr = PetscSectionSetNumFields(coordSection, 1);CHKERRQ(ierr);
  ierr = PetscSectionSetFieldComponents(coordSection, 0, dimEmbed);CHKERRQ(ierr);
  ierr = PetscSectionSetChart(coordSection, numCells, numCells + numVertices);CHKERRQ(ierr);
  for (v = numCells; v < numCells+numVertices; ++v) {
    ierr = PetscSectionSetDof(coordSection, v, dimEmbed);CHKERRQ(ierr);
    ierr = PetscSectionSetFieldDof(coordSection, v, 0, dimEmbed);CHKERRQ(ierr);
  }
  ierr = PetscSectionSetUp(coordSection);CHKERRQ(ierr);
  ierr = PetscSectionGetStorageSize(coordSection, &coordSize);CHKERRQ(ierr);
  ierr = VecCreate(comm, &coordinates);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) coordinates, "coordinates");CHKERRQ(ierr);
  ierr = VecSetSizes(coordinates, coordSize, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetBlockSize(coordinates, dimEmbed);CHKERRQ(ierr);
  ierr = VecSetType(coordinates, VECSTANDARD);CHKERRQ(ierr);
  ierr = VecGetArray(coordinates, &cds);CHKERRQ(ierr);
  for (v = 0; v < numVertices; ++v) {
    for (d = 0; d < dimEmbed; ++d) {
      cds[v*dimEmbed+d] = coords[v*dimEmbed+d];
    }
  }
  ierr = VecRestoreArray(coordinates, &cds);CHKERRQ(ierr);
  ierr = DMSetCoordinatesLocal(*dm, coordinates);CHKERRQ(ierr);
  ierr = VecDestroy(&coordinates);CHKERRQ(ierr);
  ierr = VecCreate(comm, n);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) *n, "normals");CHKERRQ(ierr);
  ierr = VecSetSizes(*n, coordSize, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetType(*n, VECSTANDARD);CHKERRQ(ierr);
  ierr = VecGetArray(*n, &cds);CHKERRQ(ierr);
  for (v = 0; v < numVertices; ++v) {
    for (d = 0; d < dimEmbed; ++d) {
      cds[v*dimEmbed+d] = normals[v*dimEmbed+d];
    }
  }
  ierr = VecRestoreArray(*n, &cds);CHKERRQ(ierr);
  if (!rank) {
    ierr = PetscFree2(coords, normals);CHKERRQ(ierr);
    ierr = PetscFree(cells);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "loadSrfIntoSurfacePoints"
/*@C
  loadSrfIntoSurfacePoints - Read a .srf file and return mesh data

  Input Parameters:
+ comm - The MPI Comm
- filename - The filename of the .srf file

  Output Parameters:
+ n    - The vertex normals
. w    - The vertex weights
. area - The panel areas
. totalArea - The total surface area
- dm   - The DM Plex object

  Level: developer

  Note:
  In order to plot this data in Matlab, use
$   v = coordinates; n = normals;
$   quiver3(v(:,1), v(:,2), v(:,3), n(:,1), n(:,2), n(:,3));

.seealso: DMPlexCreateBardhanFromFile(), DMPlexCreateBardhan()
@*/
PetscErrorCode loadSrfIntoSurfacePoints(MPI_Comm comm, const char filename[], Vec *n, Vec *w, Vec *area, PetscReal *totalArea, DM *dm)
{
  PetscViewer    viewer;
  PetscScalar   *a, *weight;
  char           basename[PETSC_MAX_PATH_LEN];
  PetscReal      totArea = 0.0;
  PetscInt       l, cStart, cEnd, vStart, vEnd, c;
  size_t         len;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscViewerCreate(comm, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIRead(viewer, /* dummy */basename, 1, NULL, PETSC_STRING);CHKERRQ(ierr);
  ierr = PetscViewerASCIIRead(viewer, /* dummy */basename, 1, NULL, PETSC_STRING);CHKERRQ(ierr);
  ierr = PetscStrcpy(basename, filename);CHKERRQ(ierr);
  ierr = PetscStrlen(basename, &len);CHKERRQ(ierr);
  for (l = len-1; basename[l] != '/' && l >= 0; --l) basename[l] = '\0';
  ierr = PetscViewerASCIIRead(viewer, &basename[l+1], 1, NULL, PETSC_STRING);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  /* This is causing gibberish output, figure out before uncommenting */
  //for(int i=0; i<PETSC_MAX_PATH_LEN; ++i)
  //printf("%c",basename[i]);
  //printf("\n\n");
  ierr = DMPlexCreateBardhanFromFile(comm, basename, PETSC_TRUE, n, dm);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(*dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(*dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  /* TODO Jay has a pass where he eliminates vertices of low weight, or which are too close to another vertex */
  ierr = DMSetFromOptions(*dm);CHKERRQ(ierr);
  ierr = DMViewFromOptions(*dm, NULL, "-dm_view");CHKERRQ(ierr);
  {
    PetscDS prob;
    PetscFE fe;
    PetscSpace sp;
    PetscInt dim, order;

    ierr = DMGetDimension(*dm, &dim);CHKERRQ(ierr);
    ierr = DMGetDS(*dm, &prob);CHKERRQ(ierr);
    ierr = PetscFECreateDefault(comm, dim, 1, PETSC_TRUE, "potential_", -1, &fe);CHKERRQ(ierr);
    ierr = PetscFEGetBasisSpace(fe, &sp);CHKERRQ(ierr);
    ierr = PetscSpaceGetDegree(sp, &order, NULL);CHKERRQ(ierr);
    if (order != 1) SETERRQ(comm, PETSC_ERR_ARG_WRONG, "Must have a linear discretization. Try using -potential_petscspace_degree 1");
    ierr = PetscDSSetDiscretization(prob, 0, (PetscObject) fe);CHKERRQ(ierr);
    ierr = PetscDSSetFromOptions(prob);CHKERRQ(ierr);
    ierr = PetscFEDestroy(&fe);CHKERRQ(ierr);
  }

  ierr = VecCreate(comm, area);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) *area, "panel areas");CHKERRQ(ierr);
  ierr = VecSetSizes(*area, cEnd-cStart, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetType(*area, VECSTANDARD);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(*dm, w);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) *w, "vertex weights");CHKERRQ(ierr);
  ierr = VecSet(*w, 0.0);CHKERRQ(ierr);

  ierr = VecGetArray(*w, &weight);CHKERRQ(ierr);
  ierr = VecGetArray(*area, &a);CHKERRQ(ierr);
  for (c = cStart; c < cEnd; ++c) {
    PetscReal area, centroid[3], normal[3];
    PetscInt *closure = NULL;
    PetscInt  clSize, cl, s = 0;

    ierr = DMPlexComputeCellGeometryFVM(*dm, c, &area, centroid, normal);CHKERRQ(ierr);
    a[c-cStart] = area;
    totArea += area;
    ierr = DMPlexGetTransitiveClosure(*dm, c, PETSC_TRUE, &clSize, &closure);CHKERRQ(ierr);
    for (cl = 0; cl < clSize*2; cl += 2) {
      const PetscInt v = closure[cl];

      if ((v < vStart) || (v >= vEnd)) continue;
      weight[v-vStart] += area/3.0;
      ++s;
    }
    ierr = DMPlexRestoreTransitiveClosure(*dm, c, PETSC_TRUE, &clSize, &closure);CHKERRQ(ierr);
    if (s != 3) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid connectivity in Bardhan mesh, %d vertices on face %d", s, c);
  }
  ierr = VecRestoreArray(*w, &weight);CHKERRQ(ierr);
  ierr = VecRestoreArray(*area, &a);CHKERRQ(ierr);
  if (totalArea) *totalArea = totArea;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "makeSphereSurface"
/*@C
  makeSphereSurface - Create a point representation of a sphereical surface

  Input Parameters:
+ comm - The MPI Comm
. origin - The center of the sphere
. radius - The radius of the sphere
- numpoints - The number of points to use

  Output Parameters:
+ w  - The vertex weights
- n  - The vertex normals
. solidangle - The angular coordinates (theta, phi)
- dm - The DM

  Level: developer

  Note: This uses the method from Park, Bardhan, Makowski, Roux (2009), following reference in there

.seealso: DMPlexCreateBardhanFromFile(), DMPlexCreateBardhan()
@*/
PetscErrorCode makeSphereSurface(MPI_Comm comm, PetscReal origin[], PetscReal radius, PetscInt numPoints, Vec *w, Vec *n, Vec *solidangle, DM *dm)
{
  const PetscReal h = 2.0/numPoints;
  PetscReal       t = -1.0 + 0.5*h;
  Vec             coordinates, angle;
  PetscScalar    *coords, *ang, *normals, *weights;
  PetscInt        i, d;
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;
  ierr = VecCreate(comm, &coordinates);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) coordinates, "coordinates");CHKERRQ(ierr);
  ierr = VecSetSizes(coordinates, numPoints*3, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetBlockSize(coordinates, 3);CHKERRQ(ierr);
  ierr = VecSetType(coordinates, VECSTANDARD);CHKERRQ(ierr);
  ierr = VecDuplicate(coordinates, n);CHKERRQ(ierr);
  ierr = VecCreate(comm, &angle);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) angle, "solid angle");CHKERRQ(ierr);
  ierr = VecSetSizes(angle, numPoints*2, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetBlockSize(angle, 2);CHKERRQ(ierr);
  ierr = VecSetType(angle, VECSTANDARD);CHKERRQ(ierr);
  ierr = VecCreate(comm, w);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) *w, "vertex weights");CHKERRQ(ierr);
  ierr = VecSetSizes(*w, numPoints, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetType(*w, VECSTANDARD);CHKERRQ(ierr);
  ierr = VecGetArray(coordinates, &coords);CHKERRQ(ierr);
  ierr = VecGetArray(*n, &normals);CHKERRQ(ierr);
  ierr = VecGetArray(angle, &ang);CHKERRQ(ierr);
  ierr = VecGetArray(*w, &weights);CHKERRQ(ierr);
  for (i = 0; i < numPoints; ++i, t += h) {
    const PetscReal theta = PetscAcosReal(t);
    const PetscReal phi   = PetscSqrtReal(PETSC_PI*numPoints) * PetscAsinReal(t);

    ang[i*2+0]     = theta;
    ang[i*2+1]     = phi;
    weights[i]     = 4.0 * PETSC_PI * PetscSqr(radius) / numPoints;
    normals[i*3+0] = PetscSinReal(theta)*PetscCosReal(phi);
    normals[i*3+1] = PetscSinReal(theta)*PetscSinReal(phi);
    normals[i*3+2] = PetscCosReal(theta);
    for (d = 0; d < 3; ++d) coords[i*3+d] = radius * normals[i*3+d] + origin[d];
  }
  ierr = VecRestoreArray(coordinates, &coords);CHKERRQ(ierr);
  ierr = VecRestoreArray(*n, &normals);CHKERRQ(ierr);
  ierr = VecRestoreArray(angle, &ang);CHKERRQ(ierr);
  ierr = VecRestoreArray(*w, &weights);CHKERRQ(ierr);
  ierr = DMCreate(comm, dm);CHKERRQ(ierr);
  ierr = DMSetType(*dm, DMSHELL);CHKERRQ(ierr);
  ierr = DMSetCoordinatesLocal(*dm, coordinates);CHKERRQ(ierr);
  ierr = VecDestroy(&coordinates);CHKERRQ(ierr);
  if (solidangle) *solidangle = angle;
  else            {ierr = VecDestroy(&angle);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscSurfaceCreateMSP"
/*@C
  PetscSurfaceCreateMSP - Create a point representation of a sphereical surface from MSP output

  Input Parameters:
+ comm - The MPI Comm
. origin - The center of the sphere
. radius - The radius of the sphere
- numpoints - The number of points to use

  Output Parameters:
+ w  - The vertex weights
- n  - The vertex normals
. solidangle - The angular coordinates (theta, phi)
- dm - The DM

  Level: developer

  Note: This uses the method from Park, Bardhan, Makowski, Roux (2009), following reference in there

.seealso: DMPlexCreateBardhanFromFile(), DMPlexCreateBardhan()
@*/
PetscErrorCode PetscSurfaceCreateMSP(MPI_Comm comm, const char pntFilename[], PetscSurface *s)
{
  PetscViewer     viewer;
  Vec             coordinates;
  PetscScalar    *c, *n, *w;
  PetscReal      *coords = NULL, *normals = NULL, *weights = NULL;
  PetscInt        Np = 0, p, d;
  PetscMPIInt     rank;
  PetscBool       flg;
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;
  ierr = PetscMemzero(s, sizeof(*s));CHKERRQ(ierr);
  ierr = PetscTestFile(pntFilename, 'r', &flg);CHKERRQ(ierr);
  if (!flg) PetscFunctionReturn(0);
  ierr = PetscViewerCreate(comm, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, pntFilename);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  if (!rank) {
    PetscInt maxSize = 1000, cnt = 3;

    ierr = PetscMalloc3(maxSize*3, &coords, maxSize*3, &normals, maxSize, &weights);CHKERRQ(ierr);
    while (cnt == 3) {
      PetscInt dummy[4];

      if (Np >= maxSize) {
        PetscReal *tmpcoords, *tmpnormals, *tmpweights;

        ierr = PetscMalloc3(maxSize*2*3, &tmpcoords, maxSize*2*3, &tmpnormals, maxSize*2, &tmpweights);CHKERRQ(ierr);
        ierr = PetscMemcpy(tmpcoords,  coords,  maxSize*3 * sizeof(PetscReal));CHKERRQ(ierr);
        ierr = PetscMemcpy(tmpnormals, normals, maxSize*3 * sizeof(PetscReal));CHKERRQ(ierr);
        ierr = PetscMemcpy(tmpweights, weights, maxSize   * sizeof(PetscReal));CHKERRQ(ierr);
        ierr = PetscFree3(coords, normals, weights);CHKERRQ(ierr);
        coords  = tmpcoords;
        normals = tmpnormals;
        weights = tmpweights;
        maxSize *= 2;
      }
      ierr = PetscViewerRead(viewer, dummy,          4, &cnt, PETSC_INT);CHKERRQ(ierr);
      ierr = PetscViewerRead(viewer, &coords[Np*3],  3, &cnt, PETSC_DOUBLE);CHKERRQ(ierr);
      ierr = PetscViewerRead(viewer, &weights[Np],   1, &cnt, PETSC_DOUBLE);CHKERRQ(ierr);
      ierr = PetscViewerRead(viewer, &normals[Np*3], 3, &cnt, PETSC_DOUBLE);CHKERRQ(ierr);
      if (cnt == 3) ++Np;
    }
  }
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  /* Create coordinates */
  ierr = VecCreate(comm, &coordinates);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) coordinates, "coordinates");CHKERRQ(ierr);
  ierr = VecSetSizes(coordinates, Np*3, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetBlockSize(coordinates, 3);CHKERRQ(ierr);
  ierr = VecSetType(coordinates, VECSTANDARD);CHKERRQ(ierr);
  /* Create normals */
  ierr = VecDuplicate(coordinates, &s->normals);CHKERRQ(ierr);
  /* Create weights */
  ierr = VecCreate(comm, &s->weights);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) s->weights, "vertex weights");CHKERRQ(ierr);
  ierr = VecSetSizes(s->weights, Np, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetType(s->weights, VECSTANDARD);CHKERRQ(ierr);
  /* Set values */
  ierr = VecGetArray(coordinates, &c);CHKERRQ(ierr);
  ierr = VecGetArray(s->normals, &n);CHKERRQ(ierr);
  ierr = VecGetArray(s->weights, &w);CHKERRQ(ierr);
  for (p = 0; p < Np; ++p) {
    w[p]     = weights[p];
    for (d = 0; d < 3; ++d) {
      c[p*3+d] = coords[p*3+d];
      n[p*3+d] = normals[p*3+d];
    }
  }
  ierr = VecRestoreArray(coordinates, &c);CHKERRQ(ierr);
  ierr = VecRestoreArray(s->normals, &n);CHKERRQ(ierr);
  ierr = VecRestoreArray(s->weights, &w);CHKERRQ(ierr);
  ierr = PetscFree3(coords, normals, weights);CHKERRQ(ierr);
  /* Create DM */
  ierr = DMCreate(comm, &s->dm);CHKERRQ(ierr);
  ierr = DMSetType(s->dm, DMSHELL);CHKERRQ(ierr);
  ierr = DMSetCoordinatesLocal(s->dm, coordinates);CHKERRQ(ierr);
  ierr = VecDestroy(&coordinates);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscSurfaceDestroy"
/*@
  PetscSurfaceDestroy - Destroy the PetscSurface

  Input Parameters:
. s - The PetscSurface

  Level: beginner

.seealso: DMPlexCreateBardhanFromFile(), DMPlexCreateBardhan()
@*/
PetscErrorCode PetscSurfaceDestroy(PetscSurface *s)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = DMDestroy(&s->dm);CHKERRQ(ierr);
  ierr = VecDestroy(&s->weights);CHKERRQ(ierr);
  ierr = VecDestroy(&s->normals);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
