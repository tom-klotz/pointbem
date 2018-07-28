#if !defined(__MOLECULE_H)
#define __MOLECULE_H

typedef struct {
  Vec q;   /* Charge values */
  Vec xyz; /* Charge coordinates, always 3D */
  Vec R;   /* Charge radii */
} PQRData;


#endif
