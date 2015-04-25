#include "meshmaker.h"
 
/* Constructors and Destructors */

Panel Panel_allocate() {
   return (Panel)calloc(1, sizeof(struct _Panel));
}

void Panel_free(Panel panel) {
/*
   switch (panel->type) {
      case FLAT_PANEL: FlatPanel_free(panel->realpanel); break;
      case GST_PANEL: GST_free(panel->realpanel); break;
      case TOR_PANEL: TOR_free(panel->realpanel); break;
   }
*/
#ifdef ENABLE_GALERKIN
  unsigned int i;
  for (i = 0; i < panel->numgalerkinquadpoints; i++) {
	 Vector3D_free(panel->galerkinquadpoints[i]);
	 Vector3D_free(panel->galerkinquadnormals[i]);
  }
  free(panel->galerkinquadnormals);
  free(panel->galerkinquadpoints);
  Vector_free(panel->galerkinquadweights);
#endif
   free(panel);
}

/* Operators */
                                                                                
void Panel_FlatPanel(Panel panel, FlatPanel flatpanel) {
   panel->centroid = flatpanel->centroid;
   panel->normal = flatpanel->panelaxis[2];

   panel->numdirectquadpoints = flatpanel->numdirectquadpoints;
   panel->directquadpoints = flatpanel->directquadpoints;
   panel->directquadnormals = flatpanel->directquadnormals;
   panel->directquadweights = flatpanel->directquadweights;

#ifdef ENABLE_GALERKIN
/* 	panel->numgalerkinquadpoints = galerkinOrder; */
/* 	panel->galerkinquadpoints = (Vector3D *)calloc(panel->numgalerkinquadpoints, sizeof(Vector3D)); */
/* 	panel->galerkinquadnormals = (Vector3D *)calloc(panel->numgalerkinquadpoints, sizeof(Vector3D)); */
/* 	panel->galerkinquadweights = Vector_allocate(panel->numgalerkinquadpoints); */
/* 	unsigned int i; */
/* 	for (i = 0; i < panel->numgalerkinquadpoints; i++) { */
/* 	  panel->galerkinquadpoints[i] = Vector3D_allocate(); */
/* 	  panel->galerkinquadnormals[i] = Vector3D_allocate(); */
/* 	} */
	FlatPanel_getquadrature(flatpanel, GALERKIN_ORDER, &(panel->numgalerkinquadpoints),
									&(panel->galerkinquadpoints), &(panel->galerkinquadnormals), &(panel->galerkinquadweights));
#endif

   Vector3D point = Vector3D_allocate();
   point->x = 1000000.0;
   //panel->area = Integration(point, panel, CONSTANT_KERNEL, NULL, SINGLE_LAYER_INT);
   panel->area = flatpanel->area;
   panel->type = FLAT_PANEL;
   panel->realpanel = flatpanel;

   Vector3D_free(point);
}

void Panel_GST(Panel panel, GST gst, unsigned int usecom) {
   panel->normal = gst->normalAtCentroid;
   panel->type = GST_PANEL;
   panel->realpanel = gst;
   panel->maxedgelength = gst->flatpanel->maxEdgelength;

   Vector3D point = Vector3D_allocate();
   point->x = 1000000.0;

   if (usecom) {
      panel->centroid = Vector3D_allocate();
      panel->centroid->x = Integration(point, panel, X_KERNEL, NULL, SINGLE_LAYER_INT) / panel->area;
      panel->centroid->y = Integration(point, panel, Y_KERNEL, NULL, SINGLE_LAYER_INT) / panel->area;
      panel->centroid->z = Integration(point, panel, Z_KERNEL, NULL, SINGLE_LAYER_INT) / panel->area;
      panel->centroid->x -= 1.0;
      panel->centroid->y -= 1.0;
      panel->centroid->z -= 1.0;
      Vector3D_sub(panel->centroid, panel->centroid, gst->Tvec);
      Vector3D_transformVecs_inverse(panel->centroid, gst->Tmat[0], gst->Tmat[1], gst->Tmat[2], panel->centroid);
   }
   else {
      panel->centroid = gst->centroid;
   }

   panel->numdirectquadpoints = gst->numdirectquadpoints;
   panel->directquadpoints = gst->directquadpoints;
   panel->directquadnormals = gst->directquadnormals;
   panel->directquadweights = gst->directquadweights;

   panel->area = Integration(point, panel, CONSTANT_KERNEL, NULL, SINGLE_LAYER_INT);

   Vector3D_free(point);
}

void Panel_TOR(Panel panel, TOR tor, unsigned int usecom) {
   panel->normal = tor->normalAtCentroid;
   panel->type = TOR_PANEL;
   panel->realpanel = tor;
   panel->maxedgelength = tor->maxEdgelength;

   Vector3D point = Vector3D_allocate();
   point->x = 1000000.0;

   if (usecom) {
      panel->centroid = Vector3D_allocate();
      panel->centroid->x = Integration(point, panel, X_KERNEL, NULL, SINGLE_LAYER_INT) / panel->area;
      panel->centroid->y = Integration(point, panel, Y_KERNEL, NULL, SINGLE_LAYER_INT) / panel->area;
      panel->centroid->z = Integration(point, panel, Z_KERNEL, NULL, SINGLE_LAYER_INT) / panel->area;
      panel->centroid->x -= 1.0;
      panel->centroid->y -= 1.0;
      panel->centroid->z -= 1.0;
      Vector3D_sub(panel->centroid, panel->centroid, tor->Tvec);
      Vector3D_transformVecs_inverse(panel->centroid, tor->Tmat[0], tor->Tmat[1], tor->Tmat[2], panel->centroid);
   }
   else {
      panel->centroid = tor->centroid;
   }

   panel->numdirectquadpoints = tor->numdirectquadpoints;
   panel->directquadpoints = tor->directquadpoints;
   panel->directquadnormals = tor->directquadnormals;
   panel->directquadweights = tor->directquadweights;

   panel->area = Integration(point, panel, CONSTANT_KERNEL, NULL, SINGLE_LAYER_INT);

   Vector3D_free(point);
}

void Panel_Vector3D(Panel panel, Vector3D v3d) {
   panel->centroid = v3d;
   panel->normal = NULL;
   panel->type = POINT_PANEL;
   panel->realpanel = v3d;
   panel->area = 0.0;
}

real Panel_memory(Panel panel, BEMLayerType layertype) {
   switch (panel->type) {
   case FLAT_PANEL: return FlatPanel_memory(4, layertype);
   case GST_PANEL: return 0;
   case TOR_PANEL: return 0;
   default: return 0;
   }
}

real Panel_quadrature(Vector3D point, Panel panel, BEMKernelType kerneltype, void* parameters, BEMLayerType layertype) {
   unsigned int i;
   real integral = 0.0;

   if (layertype == SINGLE_LAYER_INT) {
      for (i = 0; i < panel->numdirectquadpoints; i++)
         integral += panel->directquadweights[i] * GreensFunction(point, panel->directquadpoints[i], kerneltype, parameters);
   }
   else if (layertype == DOUBLE_LAYER_INT) {
      for (i = 0; i < panel->numdirectquadpoints; i++)
         integral -= panel->directquadweights[i] * GreensFunction_deriv(point, panel->directquadpoints[i], kerneltype, parameters, panel->directquadnormals[i]);
   }

   return integral;
}
