#include "meshmaker.h"

unsigned int quadcount = 0;
unsigned int allcount = 0;

inline real Integration_single_XGB(Vector3D point, Panel panel, BEMKernelType kerneltype, void *parameters, BEMLayerType layertype) {

  Vector3D ultDest = Vector3D_allocate();
  Vector3D rvector_srcToQuad = Vector3D_allocate();
  Vector3D rvector_quadToDest = Vector3D_allocate();
  FlatPanel fp = (FlatPanel)panel->realpanel;
  Vector3D normal = fp->panelaxis[2];

  ultDest->x = ((real*)parameters)[0];
  ultDest->y = ((real*)parameters)[1];
  ultDest->z = ((real*)parameters)[2];

  Vector3D_sub(rvector_srcToQuad, point, panel->centroid);
  real r_srcQuad = Vector3D_length(rvector_srcToQuad);
  Vector3D_sub(rvector_quadToDest, ultDest, panel->centroid);
  real r_quadUlt = Vector3D_length(rvector_quadToDest);

  real g = -panel->area * Vector3D_dot(rvector_srcToQuad, normal) / (pow(r_srcQuad, 3.0) * r_quadUlt);
  
  Vector3D_free(rvector_quadToDest);
  Vector3D_free(rvector_srcToQuad);
  Vector3D_free(ultDest);

  return g;
}

real Integration(Vector3D point, Panel panel, BEMKernelType kerneltype, void* parameters, BEMLayerType layertype) {
   real slp = 0.0, dlp = 0.0;

   allcount++;

   //return Panel_quadrature(point, panel, kerneltype, parameters, layertype);

/* 	if ((panel->type == FLAT_PANEL) && (kerneltype == XGB_KERNEL)) { */
/* 	  real dist = Vector3D_distance(panel->centroid, point); */
/* 	  if (dist > 10.0 * panel->maxedgelength) */
/* 		 return Integration_single_XGB(point, panel, kerneltype, parameters, layertype); */
/* 	} */
	
   if ((panel->type == GST_PANEL) || (panel->type == TOR_PANEL)) {
      real dist = Vector3D_distance(panel->centroid, point);
      if ((layertype == SINGLE_LAYER_INT) && (dist > 4.0 * panel->maxedgelength)) {
        quadcount++;
        return Panel_quadrature(point, panel, kerneltype, parameters, layertype);
      }
      else if ((layertype == DOUBLE_LAYER_INT) && (dist > 2.0 * panel->maxedgelength)) {
        quadcount++;
        return Panel_quadrature(point, panel, kerneltype, parameters, layertype);
      }
   }

   switch (panel->type) {
   case FLAT_PANEL: {
      switch (kerneltype) {
      case CONSTANT_KERNEL: slp = panel->area; break;
      case POISSON_KERNEL: FlatIntegration_oneoverr(point, (FlatPanel)(panel->realpanel), parameters, &slp, &dlp); break;
      case HELMHOLTZ_KERNEL: FlatIntegration_ekroverr_numerical(point, (FlatPanel)(panel->realpanel), parameters, &slp, &dlp); break;
      //case HELMHOLTZ_KERNEL: FlatIntegration_ekroverr_desingularized(point, (FlatPanel)(panel->realpanel), parameters, &slp, &dlp); break;
      /* These kernels are not singular, or singular but always 
         evaluated far from the panel.  Always use quadrature */
      case LJ_KERNEL: return FlatIntegration_maltquad(point, (FlatPanel)(panel->realpanel), kerneltype, parameters, layertype);
      case LJ12_KERNEL: return FlatIntegration_maltquad(point, (FlatPanel)(panel->realpanel), kerneltype, parameters, layertype);
      case LJ6_KERNEL: return FlatIntegration_maltquad(point, (FlatPanel)(panel->realpanel), kerneltype, parameters, layertype);
      case MONOMIAL_KERNEL: return FlatIntegration_maltquad(point, (FlatPanel)(panel->realpanel), kerneltype, parameters, layertype);
      case INVERSEPOWER_KERNEL: return FlatIntegration_maltquad(point, (FlatPanel)(panel->realpanel), kerneltype, parameters, layertype);
      case GHOSH_KERNEL: return FlatIntegration_maltquad(point, (FlatPanel)(panel->realpanel), kerneltype, parameters, layertype);
      case GRYCUK_KERNEL: return FlatIntegration_maltquad(point, (FlatPanel)(panel->realpanel), kerneltype, parameters, layertype);
      case LESYNG_KERNEL: return FlatIntegration_maltquad(point, (FlatPanel)(panel->realpanel), kerneltype, parameters, layertype);
      case XGB_KERNEL: return FlatIntegration_maltquad(point, (FlatPanel)(panel->realpanel), kerneltype, parameters, layertype);
      default: error("Invalid kernel type %d", kerneltype);
      }
      break;
   }
   case GST_PANEL: {
      switch (kerneltype) {
      case CONSTANT_KERNEL: Integration_general_GST_single(point, (GST)(panel->realpanel), CONSTANT_KERNEL, parameters, &slp); break;
      case X_KERNEL: Integration_general_GST_single(point, (GST)(panel->realpanel), X_KERNEL, parameters, &slp); break;
      case Y_KERNEL: Integration_general_GST_single(point, (GST)(panel->realpanel), Y_KERNEL, parameters, &slp); break;
      case Z_KERNEL: Integration_general_GST_single(point, (GST)(panel->realpanel), Z_KERNEL, parameters, &slp); break;
      case LJ12_KERNEL: Integration_general_GST_double(point, (GST)(panel->realpanel), LJ12_KERNEL, parameters, &dlp); break;
      case LJ6_KERNEL: Integration_general_GST_double(point, (GST)(panel->realpanel), LJ6_KERNEL, parameters, &dlp); break;
         //      case BCA_KERNEL: Integration_general_GST_double(point, (GST)(panel->realpanel), BCA_KERNEL, parameters, &dlp); break;
      case GRYCUK_KERNEL: dlp = Panel_quadrature(point, panel, LJ6_KERNEL, parameters, DOUBLE_LAYER_INT); break;
      case LESYNG_KERNEL: dlp = Panel_quadrature(point, panel, LESYNG_KERNEL, parameters, DOUBLE_LAYER_INT); break;
         //      case GHOSH_KERNEL: Integration_general_GST_double(point, (GST)(panel->realpanel), GHOSH_KERNEL, parameters, &dlp); break;
      case GHOSH_KERNEL: dlp = Panel_quadrature(point, panel, GHOSH_KERNEL, parameters, DOUBLE_LAYER_INT); break;
		case XGB_KERNEL: dlp = Panel_quadrature(point, panel, XGB_KERNEL, parameters, DOUBLE_LAYER_INT); break;
      case MONOMIAL_KERNEL: {
         switch (layertype) {
         case SINGLE_LAYER_INT: Integration_general_GST_single(point, (GST)(panel->realpanel), MONOMIAL_KERNEL, parameters, &slp); break;
         case DOUBLE_LAYER_INT: Integration_general_GST_double(point, (GST)(panel->realpanel), MONOMIAL_KERNEL, parameters, &dlp); break;
         default: error("Invalid layer type %d", layertype);
         }
         break;
      }
      case INVERSEPOWER_KERNEL: {
         switch (layertype) {
         case SINGLE_LAYER_INT: Integration_general_GST_single(point, (GST)(panel->realpanel), INVERSEPOWER_KERNEL, parameters, &slp); break;
         case DOUBLE_LAYER_INT: Integration_general_GST_double(point, (GST)(panel->realpanel), INVERSEPOWER_KERNEL, parameters, &dlp); break;
         default: error("Invalid layer type %d", layertype);
         }
         break;
      }
      case POISSON_KERNEL: {
         switch (layertype) {
         case SINGLE_LAYER_INT: Integration_oneoverr_GST_single(point, (GST)(panel->realpanel), parameters, &slp); break;
         case DOUBLE_LAYER_INT: Integration_oneoverr_GST_double(point, (GST)(panel->realpanel), parameters, &dlp); break;
         default: error("Invalid layer type %d", layertype);
         }
         break;
      }
      case HELMHOLTZ_KERNEL: {
         switch (layertype) {
         case SINGLE_LAYER_INT: Integration_oneoverr_GST_single(point, (GST)(panel->realpanel), parameters, &slp);
           slp += Panel_quadrature(point, panel, DESINGULARIZED_HELMHOLTZ_KERNEL, parameters, SINGLE_LAYER_INT);
           break;
         case DOUBLE_LAYER_INT: Integration_oneoverr_GST_double(point, (GST)(panel->realpanel), parameters, &dlp);
           dlp += Panel_quadrature(point, panel, DESINGULARIZED_HELMHOLTZ_KERNEL, parameters, DOUBLE_LAYER_INT);
           break;
         default: error("Invalid layer type %d", layertype);
         }
         break;
      }
      default: error("Invalid kernel type %d", kerneltype);
      }
      break;
   }
   case TOR_PANEL: {
      switch (kerneltype) {
      case CONSTANT_KERNEL: Integration_general_TOR_single(point, (TOR)(panel->realpanel), CONSTANT_KERNEL, parameters, &slp); break;
      case X_KERNEL: Integration_general_TOR_single(point, (TOR)(panel->realpanel), X_KERNEL, parameters, &slp); break;
      case Y_KERNEL: Integration_general_TOR_single(point, (TOR)(panel->realpanel), Y_KERNEL, parameters, &slp); break;
      case Z_KERNEL: Integration_general_TOR_single(point, (TOR)(panel->realpanel), Z_KERNEL, parameters, &slp); break;
      case LESYNG_KERNEL: dlp = Panel_quadrature(point, panel, LESYNG_KERNEL, parameters, DOUBLE_LAYER_INT); break;
      case LJ6_KERNEL: dlp = Panel_quadrature(point, panel, LJ6_KERNEL, parameters, DOUBLE_LAYER_INT); break;
      case GRYCUK_KERNEL: dlp = Panel_quadrature(point, panel, LJ6_KERNEL, parameters, DOUBLE_LAYER_INT); break;
      case GHOSH_KERNEL: dlp = Panel_quadrature(point, panel, GHOSH_KERNEL, parameters, DOUBLE_LAYER_INT); break;
		case XGB_KERNEL: dlp = Panel_quadrature(point, panel, XGB_KERNEL, parameters, DOUBLE_LAYER_INT); break;
      case MONOMIAL_KERNEL: {
         switch (layertype) {
         case SINGLE_LAYER_INT: Integration_general_TOR_single(point, (TOR)(panel->realpanel), MONOMIAL_KERNEL, parameters, &slp); break;
         //case DOUBLE_LAYER_INT: Integration_general_TOR_double(point, (TOR)(panel->realpanel), MONOMIAL_KERNEL, parameters, &dlp); break;
         default: error("Invalid layer type %d", layertype);
         }
         break;
      }
      case INVERSEPOWER_KERNEL: {
         switch (layertype) {
         case SINGLE_LAYER_INT: Integration_general_TOR_single(point, (TOR)(panel->realpanel), INVERSEPOWER_KERNEL, parameters, &slp); break;
         //case DOUBLE_LAYER_INT: Integration_general_TOR_double(point, (TOR)(panel->realpanel), INVERSEPOWER_KERNEL, parameters, &dlp); break;
         default: error("Invalid layer type %d", layertype);
         }
         break;
      }
      case POISSON_KERNEL: {
         switch (layertype) {
         case SINGLE_LAYER_INT: Integration_oneoverr_TOR_single(point, (TOR)(panel->realpanel), parameters, &slp); break;
         case DOUBLE_LAYER_INT: Integration_oneoverr_TOR_double(point, (TOR)(panel->realpanel), parameters, &dlp); break;
         default: error("Invalid layer type %d", layertype);
         }
         break;
      }
      case HELMHOLTZ_KERNEL: {
         switch (layertype) {
         case SINGLE_LAYER_INT: Integration_oneoverr_TOR_single(point, (TOR)(panel->realpanel), parameters, &slp);
           slp += Panel_quadrature(point, panel, DESINGULARIZED_HELMHOLTZ_KERNEL, parameters, SINGLE_LAYER_INT);
           break;
         case DOUBLE_LAYER_INT: Integration_oneoverr_TOR_double(point, (TOR)(panel->realpanel), parameters, &dlp);
           dlp += Panel_quadrature(point, panel, DESINGULARIZED_HELMHOLTZ_KERNEL, parameters, DOUBLE_LAYER_INT);
           break;
         default: error("Invalid layer type %d", layertype);
         }
         break;
      }
      default: error("Invalid kernel type %d", kerneltype);
      }
      break;
   }
   case POINT_PANEL: {
      slp = GreensFunction(point, panel->centroid, kerneltype, parameters); break;
   }
   default: error("Invalid panel type %d", panel->type);
   }

   switch (layertype) {
      case SINGLE_LAYER_INT: return slp;
      case DOUBLE_LAYER_INT: return dlp;
      default: return 0.0;
   }
}
