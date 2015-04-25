#include "meshmaker.h"

/* Constructors and Destructors */

VertFace VertFace_allocate() {
   return (VertFace)calloc(1, sizeof(_VertFace));
}

void VertFace_free(VertFace vf) {
   unsigned int i;
   
   for (i = 0; i < vf->numvertices; i++)
      Vector3D_free(vf->vertices[i]);

   free(vf->vertices);
   free(vf->facesx);
   free(vf->facesy);
   free(vf->facesz);
   free(vf->facesgenus);

   free(vf);
}

/* Operations */

void VertFace_readvert(VertFace vf, FILE* file) {
   char line[128];
   unsigned int i;
        
   vf->numvertices = 0;
    
   while (fgets(line, 128, file))
      vf->numvertices++;
       
   vf->vertices = (Vector3D*)calloc(vf->numvertices, sizeof(Vector3D));
    
   rewind(file);
 
   for (i = 0; i < vf->numvertices; i++) {
      vf->vertices[i] = Vector3D_allocate();
   
      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf", &vf->vertices[i]->x, &vf->vertices[i]->y, &vf->vertices[i]->z);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f", &vf->vertices[i]->x, &vf->vertices[i]->y, &vf->vertices[i]->z);
#endif
#endif
   }
}
 
void VertFace_readface(VertFace vf, FILE* file) {
   char line[128];
   unsigned int i;
        
   vf->numfaces = 0;
    
   while (fgets(line, 128, file))
      vf->numfaces++;
       
   vf->facesx = (unsigned int*)calloc(vf->numfaces, sizeof(unsigned int));
   vf->facesy = (unsigned int*)calloc(vf->numfaces, sizeof(unsigned int));
   vf->facesz = (unsigned int*)calloc(vf->numfaces, sizeof(unsigned int));
   vf->facesgenus = (unsigned int*)calloc(vf->numfaces, sizeof(unsigned int));
    
   rewind(file);
 
   for (i = 0; i < vf->numfaces; i++) {
      fgets(line, 128, file);
      sscanf(line, "%u %u %u %u", &vf->facesx[i], &vf->facesy[i], &vf->facesz[i], &vf->facesgenus[i]);
      vf->facesx[i]--;
      vf->facesy[i]--;
      vf->facesz[i]--;
   }
}

void VertFace_readface_flip(VertFace vf, FILE* file) {
   char line[128];
   unsigned int i;
        
   vf->numfaces = 0;
    
   while (fgets(line, 128, file))
      vf->numfaces++;
       
   vf->facesx = (unsigned int*)calloc(vf->numfaces, sizeof(unsigned int));
   vf->facesy = (unsigned int*)calloc(vf->numfaces, sizeof(unsigned int));
   vf->facesz = (unsigned int*)calloc(vf->numfaces, sizeof(unsigned int));
   vf->facesgenus = (unsigned int*)calloc(vf->numfaces, sizeof(unsigned int));
    
   rewind(file);
 
   for (i = 0; i < vf->numfaces; i++) {
      fgets(line, 128, file);
      sscanf(line, "%u %u %u %u", &vf->facesx[i], &vf->facesz[i], &vf->facesy[i], &vf->facesgenus[i]);
      vf->facesx[i]--;
      vf->facesy[i]--;
      vf->facesz[i]--;
   }
}

void VertFace_fix(VertFace vf, unsigned int accessible) {
   /* This function gets rid of invalid faces */

   unsigned int goodfaces = 0;
   unsigned int* newfacesx;
   unsigned int* newfacesy;
   unsigned int* newfacesz;
   unsigned int* newfacesgenus;

   unsigned int f;

   for (f = 0; f < vf->numfaces; f++) {
      real d1 = Vector3D_distance(vf->vertices[vf->facesx[f]], vf->vertices[vf->facesy[f]]);
      real d2 = Vector3D_distance(vf->vertices[vf->facesy[f]], vf->vertices[vf->facesz[f]]);
      real d3 = Vector3D_distance(vf->vertices[vf->facesz[f]], vf->vertices[vf->facesx[f]]);
      real s = 0.5 * (d1 + d2 + d3);
      real area = sqrt(s * (s - d1) * (s - d2) * (s - d3));
      if (area > 0.0) {
#ifndef SCATTER
         if (accessible) {  /* Accessible surfaces should only have caps */
            if (vf->facesgenus[f] == 1)  /* This is a cap face */
               goodfaces++;
         }
         else
            goodfaces++;
#else
			goodfaces++;
#endif
      }
   }

   newfacesx = (unsigned int*)calloc(goodfaces, sizeof(unsigned int));
   newfacesy = (unsigned int*)calloc(goodfaces, sizeof(unsigned int));
   newfacesz = (unsigned int*)calloc(goodfaces, sizeof(unsigned int));
   newfacesgenus = (unsigned int*)calloc(goodfaces, sizeof(unsigned int));

   unsigned int goodcount = 0;

   for (f = 0; f < vf->numfaces; f++) {
      real d1 = Vector3D_distance(vf->vertices[vf->facesx[f]], vf->vertices[vf->facesy[f]]);
      real d2 = Vector3D_distance(vf->vertices[vf->facesy[f]], vf->vertices[vf->facesz[f]]);
      real d3 = Vector3D_distance(vf->vertices[vf->facesz[f]], vf->vertices[vf->facesx[f]]);
      real s = 0.5 * (d1 + d2 + d3);
      real area = sqrt(s * (s - d1) * (s - d2) * (s - d3));
      if (area > 0.0) {
#ifndef SCATTER
		  if (accessible) {  /* Accessible surfaces should only have caps */
            if (vf->facesgenus[f] == 1) { /* This is a cap face */
               newfacesx[goodcount] = vf->facesx[f];
               newfacesy[goodcount] = vf->facesy[f];
               newfacesz[goodcount] = vf->facesz[f];
               newfacesgenus[goodcount] = vf->facesgenus[f];
               goodcount++;
            }
         }
         else {
            newfacesx[goodcount] = vf->facesx[f];
            newfacesy[goodcount] = vf->facesy[f];
            newfacesz[goodcount] = vf->facesz[f];
            newfacesgenus[goodcount] = vf->facesgenus[f];
            goodcount++;
         }
#else
		  newfacesx[goodcount] = vf->facesx[f];
		  newfacesy[goodcount] = vf->facesy[f];
		  newfacesz[goodcount] = vf->facesz[f];
		  newfacesgenus[goodcount] = vf->facesgenus[f];
		  goodcount++;
#endif
      } else {
		  printf("face %d marked as bad!\n", f);
		}
   }

   free(vf->facesx);
   free(vf->facesy);
   free(vf->facesz);
   free(vf->facesgenus);

   vf->facesx = newfacesx;
   vf->facesy = newfacesy;
   vf->facesz = newfacesz;
   vf->facesgenus = newfacesgenus;

   printf("Removed %u panels with zero area or invalid genus\n", vf->numfaces - goodfaces);

   vf->numfaces = goodfaces;
}

void VertFace_getpanels(VertFace vf, Panel* panels) {
   unsigned int i;
    
   for (i = 0; i < vf->numfaces; i++) {
      FlatPanel flatpanel = FlatPanel_allocate(vf->vertices[vf->facesx[i]],
                                               vf->vertices[vf->facesy[i]],
                                               vf->vertices[vf->facesz[i]]);
      panels[i] = Panel_allocate();
      Panel_FlatPanel(panels[i], flatpanel);
   }
}
