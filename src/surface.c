#include <petsc.h>
#include <constants.h>
#include <surface.h>

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
    PetscInt maxSize = 1000, cnt = 3, tot = 0;

    ierr = PetscMalloc2(maxSize, &coords, maxSize, &normals);CHKERRQ(ierr);
    while (cnt == 3) {
      PetscInt dummy[3];

      if (tot + cnt > maxSize) {
        PetscReal *tmpcoords, *tmpnormals;

        ierr = PetscMalloc2(maxSize*2, &tmpcoords, maxSize*2, &tmpnormals);CHKERRQ(ierr);
        ierr = PetscMemcpy(tmpcoords,  coords,  maxSize * sizeof(PetscReal));CHKERRQ(ierr);
        ierr = PetscMemcpy(tmpnormals, normals, maxSize * sizeof(PetscReal));CHKERRQ(ierr);
        coords  = tmpcoords;
        normals = tmpnormals;
        ierr = PetscFree2(coords, normals);CHKERRQ(ierr);
        maxSize *= 2;
      }
      ierr = PetscViewerRead(viewerVert, &coords[tot],  &cnt, PETSC_FLOAT);CHKERRQ(ierr);
      ierr = PetscViewerRead(viewerVert, &normals[tot], &cnt, PETSC_FLOAT);CHKERRQ(ierr);
      ierr = PetscViewerRead(viewerVert, dummy,         &cnt, PETSC_INT);CHKERRQ(ierr);
      tot += cnt;
    }
    numVertices = tot;

    cnt  = 3;
    tot  = 3;
    ierr = PetscMalloc1(maxSize, &cells);CHKERRQ(ierr);
    while (cnt == 3) {
      PetscInt dummy[3];

      if (tot + cnt > maxSize) {
        PetscInt *tmpcells;

        ierr = PetscMalloc1(maxSize*2, &tmpcells);CHKERRQ(ierr);
        ierr = PetscMemcpy(tmpcells, cells, maxSize * sizeof(PetscInt));CHKERRQ(ierr);
        cells = tmpcells;
        ierr = PetscFree(cells);CHKERRQ(ierr);
        maxSize *= 2;
      }
      ierr = PetscViewerRead(viewerFace, &cells[tot], &cnt, PETSC_INT);CHKERRQ(ierr);
      ierr = PetscViewerRead(viewerFace, dummy,       &cnt, PETSC_INT);CHKERRQ(ierr);
      tot += cnt;
    }
    numCells = tot;
  }

  /* Allocate cell-vertex mesh */
  ierr = DMCreate(comm, dm);CHKERRQ(ierr);
  ierr = DMSetType(*dm, DMPLEX);CHKERRQ(ierr);
  ierr = DMSetDimension(*dm, dim);CHKERRQ(ierr);
  ierr = DMPlexSetChart(*dm, 0, numCells + numVertices);CHKERRQ(ierr);
  if (numCells < 0) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid number of cells %d in Bardhan file", numCells);
  for (c = 0; c < numCells; ++c) {ierr = DMPlexSetConeSize(*dm, c, numCellVertices);CHKERRQ(ierr);}
  ierr = DMSetUp(*dm);CHKERRQ(ierr);
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
  ierr = VecSetType(coordinates, VECSTANDARD);CHKERRQ(ierr);
  ierr = VecGetArray(coordinates, &cds);CHKERRQ(ierr);
  for (v = 0; v < numVertices; ++v) {
    for (d = 0; d < dimEmbed; ++d) {
      cds[v*dim+d] = coords[v*dim+d];
    }
  }
  ierr = VecRestoreArray(coordinates, &cds);CHKERRQ(ierr);
  ierr = DMSetCoordinatesLocal(*dm, coordinates);CHKERRQ(ierr);
  ierr = VecDestroy(&coordinates);CHKERRQ(ierr);
  ierr = VecCreate(comm, n);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) coordinates, "normals");CHKERRQ(ierr);
  ierr = VecSetSizes(*n, coordSize, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetType(*n, VECSTANDARD);CHKERRQ(ierr);
  ierr = VecGetArray(*n, &cds);CHKERRQ(ierr);
  for (v = 0; v < numVertices; ++v) {
    for (d = 0; d < dimEmbed; ++d) {
      cds[v*dim+d] = normals[v*dim+d];
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
+ n  - The vertex normals
. w  - The vertex weights
- dm - The DM Plex object

  Level: developer

  Note:
  In order to plot this data in Matlab, use
$   v = coordinates; n = normals;
$   quiver3(v(:,1), v(:,2), v(:,3), n(:,1), n(:,2), n(:,3));

.seealso: DMPlexCreateBardhanFromFile(), DMPlexCreateBardhan()
@*/
PetscErrorCode loadSrfIntoSurfacePoints(MPI_Comm comm, const char filename[], Vec *n, Vec *w, DM *dm)
{
  PetscViewer    viewer;
  PetscScalar   *a;
  char           basename[PETSC_MAX_PATH_LEN];
  PetscInt       num = PETSC_MAX_PATH_LEN-1, cStart, cEnd, vStart, vEnd, c;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscViewerCreate(comm, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filename);CHKERRQ(ierr);
  ierr = PetscViewerRead(viewer, /* dummy */basename, &num, PETSC_STRING);CHKERRQ(ierr);
  ierr = PetscViewerRead(viewer, /* dummy */basename, &num, PETSC_STRING);CHKERRQ(ierr);
  ierr = PetscViewerRead(viewer, /* dummy */basename, &num, PETSC_STRING);CHKERRQ(ierr);
  ierr = PetscViewerRead(viewer, basename, &num, PETSC_STRING);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);

  ierr = DMPlexCreateBardhanFromFile(comm, basename, PETSC_TRUE, n, dm);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(*dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(*dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  /* TODO Jay has a pass where he eliminates vertices of low weight, or which are too close to another vertex */

  ierr = VecCreate(comm, w);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) *w, "vertex weights");CHKERRQ(ierr);
  ierr = VecSetSizes(*w, vEnd-vStart, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetType(*w, VECSTANDARD);CHKERRQ(ierr);
  ierr = VecSet(*w, 0.0);CHKERRQ(ierr);
  ierr = VecGetArray(*w, &a);CHKERRQ(ierr);
  for (c = cStart; c < cEnd; ++c) {
    PetscReal area, centroid[3], normal[3];
    PetscInt *closure = NULL;
    PetscInt  clSize, cl, s = 0;

    ierr = DMPlexComputeCellGeometryFVM(*dm, c, &area, centroid, normal);CHKERRQ(ierr);
    ierr = DMPlexGetTransitiveClosure(*dm, c, PETSC_TRUE, &clSize, &closure);CHKERRQ(ierr);
    for (cl = 0; cl < clSize*2; cl += 2) {
      const PetscInt v = closure[cl];

      if ((v < vStart) || (v >= vEnd)) continue;
      a[v] += area/3.0;
      ++s;
    }
    ierr = DMPlexRestoreTransitiveClosure(*dm, c, PETSC_TRUE, &clSize, &closure);CHKERRQ(ierr);
    if (s != 3) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid connectivity in Bardhan mesh, %d vertices on face %d", s, c);
  }
  ierr = VecRestoreArray(*w, &a);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
