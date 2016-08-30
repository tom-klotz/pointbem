#if !defined(__ELLIPSOID_H)
#define __ELLIPSOID_H

typedef struct {
  PetscReal a, b, c;
  PetscReal origin[3];
  PetscReal rotation[3];
  
} Ellipsoid;

#endif
