#include "meshmaker.h"


#define MESHSES_PATH "/Users/jbardhan/bin/frentos/meshses3.prl"
#define MESHSAS_PATH "/Users/jbardhan/bin/frentos/meshsas3.prl"
#define REMESH_PATH "/Users/jbardhan/bin/frentos/remesh.prl"

#define CONNECTED_GST_TOR_PATH "/Users/jbardhan/bin/frentos/connectedgsttor"

#define DIELECTRIC_NONE 0
#define DIELECTRIC_MSMS_FLAT 1
#define DIELECTRIC_MSMS_CURVED 2
#define DIELECTRIC_NETGEN_FLAT 3  /* does not exist */
#define DIELECTRIC_NETGEN_CURVED 4

#define SALT_NONE 0
#define SALT_MSMS_FLAT 1
#define SALT_MSMS_CURVED 2
#define SALT_NETGEN_FLAT 3
#define SALT_NETGEN_CURVED 4

//#include <config.h>

PDBentry* PDBentries = NULL;
unsigned int numPDBentries;
SIZentry* SIZentries = NULL;
unsigned int numSIZentries;

real dielectricdensity = 1.0;
real saltdensity = 1.0;
unsigned int dielectriccode = 0;
unsigned int saltcode = 0;

unsigned int numsalts = 0;
Panel** saltpanels = NULL;
unsigned int* numsaltpanels = NULL;
char** saltfilenames = NULL;

unsigned int numdielectrics = 0;
Panel** dielectricpanels = NULL;
unsigned int* numdielectricpanels = NULL;
char** dielectricfilenames = NULL;

unsigned int numdielectriccavities = 0;
Panel** dielectriccavitypanels = NULL;
unsigned int* numdielectriccavitypanels = NULL;
char** dielectriccavityfilenames = NULL;

unsigned int numsaltcavities = 0;
Panel** saltcavitypanels = NULL;
unsigned int* numsaltcavitypanels = NULL;
char** saltcavityfilenames = NULL;

unsigned int numtotalpanels = 0;
unsigned int num_GMRES_iter = 0;
char prefix[1024];

unsigned int isInsideExhaustive(Panel* inpanels, unsigned int numinpanels, Panel* outpanels, unsigned int numoutpanels) {
   real sum = 0.0;
   unsigned int i, j;

   for (i = 0; i < numoutpanels; i++)
      for (j = 0; j < numinpanels; j++)
         sum += Integration(inpanels[j]->centroid, outpanels[i], POISSON_KERNEL, NULL, DOUBLE_LAYER_INT);

   printf("EXHAUSTIVE SUM: %f\n", sum / numinpanels);

   if (fabs(sum + 4.0*M_PI*numinpanels) < (1e-1*numinpanels))
      return 1;
   else if (fabs(sum - 4.0*M_PI*numinpanels) < (1e-1*numinpanels))
      return 1;
   else if (fabs(sum) < (1e-1*numinpanels))
      return 0;
   else {
      printf("isInsideExhaustive Futz!\n");
      return 0;
   }
}

unsigned int isInside(Panel* inpanels, unsigned int numinpaanels, Panel* outpanels, unsigned int numoutpanels) {
   real sum = 0.0;
   unsigned int i;

   for (i = 0; i < numoutpanels; i++)
      sum += Integration(inpanels[0]->centroid, outpanels[i], POISSON_KERNEL, NULL, DOUBLE_LAYER_INT);

   printf("SUM: %f\n", sum);

   if (fabs(sum + 4.0*M_PI) < 1e-1)
      return 1;
   else if (fabs(sum - 4.0*M_PI) < 1e-1)
      return 1;
   else if (fabs(sum) < 1e-1)
      return 0;
   else {
      printf("isInside Futz!\n");
      return 0;
   }
}

void generateXYZR(unsigned int accessible, char* xyzrfilename) {
   unsigned int a;
   FILE* xyzrfile = NULL;

   strcpy(xyzrfilename, "/tmp/xyzrXXXXXX");
   xyzrfile = fdopen(mkstemp(xyzrfilename), "w");

   /* REMOVE INSIDE ATOMS! */

   for (a = 0; a < numPDBentries; a++) {
      if (PDBentries[a].radius > 0.0)
            fprintf(xyzrfile, "%10.5f %10.5f %10.5f %10.5f\n",
                    PDBentries[a].x, PDBentries[a].y, PDBentries[a].z,
                    PDBentries[a].radius + (accessible ? ionexclusionradius : 0.0));
      }

   fclose(xyzrfile);
}

void generateMSMSflat(unsigned int accessible, Panel*** surfpanels, unsigned int** numsurfpanels, unsigned int* numsurfs, char*** surffilenames, Panel*** cavitypanels, unsigned int** numcavitypanels, unsigned int* numcavities, char*** cavityfilenames) {
   char xyzrfilename[256], surfacefilename[1024];
   char outputfilename[] = "/tmp/outXXXXXX";
   char command[1024];

   sprintf(surfacefilename, "%s/surfXXXXXX", prefix);

   generateXYZR(accessible, xyzrfilename);

   close(mkstemp(surfacefilename));
   close(mkstemp(outputfilename));
	char *msmspath = getenv("MSMS_PATH");

   sprintf(command, "%s -no_header -all_components -probe_radius %f -density %f -if %s -of %s > %s",
           msmspath, (accessible ? 0.002 : proberadius), (accessible ? saltdensity : dielectricdensity), xyzrfilename, surfacefilename, outputfilename);
	printf("command = %s\n", command);
	printf("output file from MSMS is %s\n", outputfilename);
	//	exit(-1);
   system(command);

   FILE* outputfile = fopen(outputfilename, "r");

   if (!outputfile)
      error("Could not open MSMS output file %s: %s", outputfilename, strerror(errno));

   char line[256];

   while (fgets(line, 256, outputfile))
      if (!strncmp(line, "NUMERICAL", 9))
         break;

   fgets(line, 256, outputfile);

   while (fgets(line, 256, outputfile)) {
      if (!strncmp(line, "    Total", 9))
         break;

      unsigned int component;
      float pr, volume, area;

      sscanf(line, "%u %f %f %f", &component, &pr, &volume, &area);
      if (volume > 0.0)
         (*numsurfs)++;
      else
         (*numcavities)++;
   }

   if ((*numsurfs) == 0)
      error("No surface with positive volume was created");

   *surfpanels = (Panel**)calloc(*numsurfs, sizeof(Panel*));
   *numsurfpanels = (unsigned int*)calloc(*numsurfs, sizeof(unsigned int));
   *surffilenames = (char**)calloc(*numsurfs, sizeof(char*));

   if ((*numcavities) > 0) {
      *cavitypanels = (Panel**)calloc(*numcavities, sizeof(Panel*));
      *numcavitypanels = (unsigned int*)calloc(*numcavities, sizeof(unsigned int));
      *cavityfilenames = (char**)calloc(*numcavities, sizeof(char*));
   }

   rewind(outputfile);

   while (fgets(line, 256, outputfile))
      if (!strncmp(line, "NUMERICAL", 9))
         break;

   fgets(line, 256, outputfile);

   unsigned int currentsurf = 0, currentcavity = 0;

   while (fgets(line, 256, outputfile)) {
      if (!strncmp(line, "    Total", 9))
         break;

      unsigned int component;
      float pr, volume, area;

      sscanf(line, "%u %f %f %f", &component, &pr, &volume, &area);

      char vertfilename[1024], facefilename[1024];

      if (component == 0) {
         sprintf(vertfilename, "%s.vert", surfacefilename);
         sprintf(facefilename, "%s.face", surfacefilename);
      }
      else {
         sprintf(vertfilename, "%s_%u.vert", surfacefilename, component);
         sprintf(facefilename, "%s_%u.face", surfacefilename, component);
      }

      FILE* vertfile = fopen(vertfilename, "r");

      if (!vertfile)
         error("Could not open MSMS vert file %s: %s", vertfilename, strerror(errno));

      FILE* facefile = fopen(facefilename, "r");

      if (!facefile)
         error("Could not open MSMS face file %s: %s", facefilename, strerror(errno));

      VertFace vertface = VertFace_allocate();

      VertFace_readvert(vertface, vertfile);

      if (volume > 0.0)
         VertFace_readface_flip(vertface, facefile);
      else
         VertFace_readface(vertface, facefile);

      VertFace_fix(vertface, accessible);

      fclose(vertfile);
      fclose(facefile);

      if (volume > 0.0) {
         (*numsurfpanels)[currentsurf] = vertface->numfaces;
         (*surfpanels)[currentsurf] = (Panel*)calloc((*numsurfpanels)[currentsurf], sizeof(Panel));
         (*surffilenames)[currentsurf] = (char*)calloc(1024, sizeof(char));

         VertFace_getpanels(vertface, (*surfpanels)[currentsurf]);

         if (component == 0)
            sprintf((*surffilenames)[currentsurf], "%s", surfacefilename);
         else
            sprintf((*surffilenames)[currentsurf], "%s_%u", surfacefilename, component);

         currentsurf++;
      }
      else {
         (*numcavitypanels)[currentcavity] = vertface->numfaces;
         (*cavitypanels)[currentcavity] = (Panel*)calloc((*numcavitypanels)[currentcavity], sizeof(Panel));
         (*cavityfilenames)[currentcavity] = (char*)calloc(1024, sizeof(char));

         VertFace_getpanels(vertface, (*cavitypanels)[currentcavity]);

         if (component == 0)
            sprintf((*cavityfilenames)[currentcavity], "%s", surfacefilename);
         else
            sprintf((*cavityfilenames)[currentcavity], "%s_%u", surfacefilename, component);

         currentcavity++;
      }

      VertFace_free(vertface);
   }

   fclose(outputfile);
/*    unlink(outputfilename); */
/*    unlink(surfacefilename); */
/*    unlink(xyzrfilename); */
}

void generateSurface_salt_MSMS_flat() {
   generateMSMSflat(1, &saltpanels, &numsaltpanels, &numsalts, &saltfilenames, 
                       &saltcavitypanels, &numsaltcavitypanels, &numsaltcavities, &saltcavityfilenames);
}

void generateSurface_salt_MSMS_curved() {
   printf("MSMS can not be used for curved salt surfaces at the current time\n");
   exit(-3);
}

void generateSurface_salt_NETGEN_flat() {
   printf("NETGEN flat not implemented\n");
   exit(-3);
}

void generateSurface_salt_NETGEN_curved() {
   char xyzrfilename[256], sasoutfilename[1024];
   char meshoutfilename[] = "/tmp/meshoutXXXXXX";
   char command[1024];

   sprintf(sasoutfilename, "%s/sasXXXXXX", prefix);

   generateXYZR(1, xyzrfilename);

   close(mkstemp(meshoutfilename));
   close(mkstemp(sasoutfilename));

   sprintf(command, "%s %s | tee %s", MESHSAS_PATH, xyzrfilename, meshoutfilename);

   system(command);

   sprintf(command, "%s %u 1 `tail -1 %s` %s", REMESH_PATH, (unsigned int)saltdensity, meshoutfilename, sasoutfilename);

   system(command);

   sprintf(command, "%s %s %s.gst /dev/null", CONNECTED_GST_TOR_PATH, xyzrfilename,  sasoutfilename);

   system(command);

   unsigned int s;

   for (s = 0; ; s++) {
      char filename[1024];
      sprintf(filename, "%s_surf_%u.gst", sasoutfilename, s);
      if (!access(filename, R_OK))
         numsalts++;
      else
         break;
   }

   for (s = 0; ; s++) {
      char filename[1024];
      sprintf(filename, "%s_cav_%u.gst", sasoutfilename, s);
      if (!access(filename, R_OK))
         numsaltcavities++;
      else
         break;
   }

   if (numsalts == 0)
      error("No external surface was created");

   saltpanels = (Panel**)calloc(numsalts, sizeof(Panel*));
   numsaltpanels = (unsigned int*)calloc(numsalts, sizeof(unsigned int));
   saltfilenames = (char**)calloc(numsalts, sizeof(char*));

   if (numsaltcavities > 0) {
      saltcavitypanels = (Panel**)calloc(numsaltcavities, sizeof(Panel*));
      numsaltcavitypanels = (unsigned int*)calloc(numsaltcavities, sizeof(unsigned int));
      saltcavityfilenames = (char**)calloc(numsaltcavities, sizeof(char*));
   }

   for (s = 0; s < numsalts; s++) {
      char filename[1024];

      sprintf(filename, "%s_surf_%u.gst", sasoutfilename, s);

      FILE* gstfile = fopen(filename, "r");

      unsigned int numgstpanels;
      GST* gstpanels;

      GST_readfile(&numgstpanels, &gstpanels, gstfile, 0);

      fclose(gstfile);

      unsigned int numtotalpanels = numgstpanels;

      numsaltpanels[s] = numtotalpanels;

      saltfilenames[s] = (char*)calloc(1024, sizeof(char));

      sprintf(saltfilenames[s], "%s_surf_%u", sasoutfilename, s);

      saltpanels[s] = (Panel*)calloc(numtotalpanels, sizeof(Panel));

      unsigned int i;

      for (i = 0; i < numgstpanels; i++) {
         saltpanels[s][i] = Panel_allocate();
         Panel_GST(saltpanels[s][i], gstpanels[i], 0);
      }
   }

   for (s = 0; s < numsaltcavities; s++) {
      char filename[1024];

      sprintf(filename, "%s_cav_%u.gst", sasoutfilename, s);

      FILE* gstfile = fopen(filename, "r");

      unsigned int numgstpanels;
      GST* gstpanels;

      GST_readfile(&numgstpanels, &gstpanels, gstfile, 0);

      fclose(gstfile);

      unsigned int numtotalpanels = numgstpanels;

      numsaltcavitypanels[s] = numtotalpanels;

      saltcavityfilenames[s] = (char*)calloc(1024, sizeof(char));

      sprintf(saltcavityfilenames[s], "%s_cav_%u", sasoutfilename, s);

      saltcavitypanels[s] = (Panel*)calloc(numtotalpanels, sizeof(Panel));

      unsigned int i;

      for (i = 0; i < numgstpanels; i++) {
         saltcavitypanels[s][i] = Panel_allocate();
         Panel_GST(saltcavitypanels[s][i], gstpanels[i], 0);
      }
   }

   char filename[1024];
   sprintf(filename, "%s.gst", sasoutfilename);
   unlink(filename);
   unlink(sasoutfilename);
}

void generateSurface_dielectric_MSMS_flat() {
   generateMSMSflat(0, &dielectricpanels, &numdielectricpanels, &numdielectrics, &dielectricfilenames, 
                       &dielectriccavitypanels, &numdielectriccavitypanels, &numdielectriccavities, &dielectriccavityfilenames);
}

void generateSurface_dielectric_MSMS_curved() {
}

void generateSurface_dielectric_NETGEN_flat() {
   printf("NETGEN can not be used for flat dielectric surfaces at the current time\n");
   exit(-3);
}

void generateSurface_dielectric_NETGEN_curved() {
   char xyzrfilename[256], sesoutfilename[1024];
   char meshoutfilename[] = "/tmp/meshoutXXXXXX";
   char command[1024];

   sprintf(sesoutfilename, "%s/sesXXXXXX", prefix);

   generateXYZR(0, xyzrfilename);

   close(mkstemp(meshoutfilename));
   close(mkstemp(sesoutfilename));

   sprintf(command, "%s %s | tee %s", MESHSES_PATH, xyzrfilename, meshoutfilename);

   system(command);

   sprintf(command, "%s %u 1 `tail -1 %s` %s", REMESH_PATH, (unsigned int)dielectricdensity, meshoutfilename, sesoutfilename);

   system(command);

   sprintf(command, "%s %s %s.gst %s.tor", CONNECTED_GST_TOR_PATH, xyzrfilename, sesoutfilename, sesoutfilename);

   system(command);

   unsigned int s;

   for (s = 0; ; s++) {
      char filename[1024];
      sprintf(filename, "%s_surf_%u.gst", sesoutfilename, s);
      if (!access(filename, R_OK))
         numdielectrics++;
      else
         break;
   }

   for (s = 0; ; s++) {
      char filename[1024];
      sprintf(filename, "%s_cav_%u.gst", sesoutfilename, s);
      if (!access(filename, R_OK))
         numdielectriccavities++;
      else
         break;
   }

   if (numdielectrics == 0)
      error("No external surface was created");

   dielectricpanels = (Panel**)calloc(numdielectrics, sizeof(Panel*));
   numdielectricpanels = (unsigned int*)calloc(numdielectrics, sizeof(unsigned int));
   dielectricfilenames = (char**)calloc(numdielectrics, sizeof(char*));

   if (numdielectriccavities > 0) {
      dielectriccavitypanels = (Panel**)calloc(numdielectriccavities, sizeof(Panel*));
      numdielectriccavitypanels = (unsigned int*)calloc(numdielectriccavities, sizeof(unsigned int));
      dielectriccavityfilenames = (char**)calloc(numdielectriccavities, sizeof(char*));
   }

   for (s = 0; s < numdielectrics; s++) {
      char filename[1024];

      sprintf(filename, "%s_surf_%u.gst", sesoutfilename, s);

      FILE* gstfile = fopen(filename, "r");

      if (!gstfile)
         error("Could not open GST file %s: %s", filename, strerror(errno));

      unsigned int numgstpanels;
      GST* gstpanels;

      GST_readfile(&numgstpanels, &gstpanels, gstfile, 0);

      fclose(gstfile);

      sprintf(filename, "%s_surf_%u.tor", sesoutfilename, s);

      FILE* torfile = fopen(filename, "r");

      if (!torfile)
         error("Could not open TOR file %s: %s", filename, strerror(errno));

      unsigned int numtorpanels;
      TOR* torpanels;

      TOR_readfile(&numtorpanels, &torpanels, torfile, 0);

      fclose(torfile);

      unsigned int numtotalpanels = numgstpanels + numtorpanels;

      numdielectricpanels[s] = numtotalpanels;

      dielectricfilenames[s] = (char*)calloc(1024, sizeof(char));

      sprintf(dielectricfilenames[s], "%s_surf_%u", sesoutfilename, s);

      dielectricpanels[s] = (Panel*)calloc(numtotalpanels, sizeof(Panel));

      unsigned int i;

      for (i = 0; i < numgstpanels; i++) {
         dielectricpanels[s][i] = Panel_allocate();
         Panel_GST(dielectricpanels[s][i], gstpanels[i], 0);
      }

      for (i = 0; i < numtorpanels; i++) {
         dielectricpanels[s][i+numgstpanels] = Panel_allocate();
         Panel_TOR(dielectricpanels[s][i+numgstpanels], torpanels[i], 0);
      }
   }

   for (s = 0; s < numdielectriccavities; s++) {
      char filename[1024];

      sprintf(filename, "%s_cav_%u.gst", sesoutfilename, s);

      FILE* gstfile = fopen(filename, "r");

      if (!gstfile)
         error("Could not open GST file %s: %s", filename, strerror(errno));

      unsigned int numgstpanels;
      GST* gstpanels;

      GST_readfile(&numgstpanels, &gstpanels, gstfile, 0);

      fclose(gstfile);

      sprintf(filename, "%s_cav_%u.tor", sesoutfilename, s);

      FILE* torfile = fopen(filename, "r");

      if (!torfile)
         error("Could not open TOR file %s: %s", filename, strerror(errno));

      unsigned int numtorpanels;
      TOR* torpanels;

      TOR_readfile(&numtorpanels, &torpanels, torfile, 0);

      fclose(torfile);

      unsigned int numtotalpanels = numgstpanels + numtorpanels;

      numdielectriccavitypanels[s] = numtotalpanels;

      dielectriccavityfilenames[s] = (char*)calloc(1024, sizeof(char));

      sprintf(dielectriccavityfilenames[s], "%s_cav_%u", sesoutfilename, s);

      dielectriccavitypanels[s] = (Panel*)calloc(numtotalpanels, sizeof(Panel));

      unsigned int i;

      for (i = 0; i < numgstpanels; i++) {
         dielectriccavitypanels[s][i] = Panel_allocate();
         Panel_GST(dielectriccavitypanels[s][i], gstpanels[i], 0);
      }

      for (i = 0; i < numtorpanels; i++) {
         dielectriccavitypanels[s][i+numgstpanels] = Panel_allocate();
         Panel_TOR(dielectriccavitypanels[s][i+numgstpanels], torpanels[i], 0);
      }
   }

   char filename[1024];
   sprintf(filename, "%s.gst", sesoutfilename);
   unlink(filename);
   sprintf(filename, "%s.tor", sesoutfilename);
   unlink(filename);
   unlink(sesoutfilename);
}

void writeSRF(char* filename) {
   unsigned int s, d, dc, sc;
   FILE* srf = fopen(filename, "w");

   if ((saltcode == SALT_MSMS_FLAT) || (saltcode == SALT_NETGEN_FLAT))
      fprintf(srf, "f\n");
   else
      fprintf(srf, "c\n");

   if ((dielectriccode == DIELECTRIC_MSMS_FLAT) || (dielectriccode == DIELECTRIC_NETGEN_FLAT))
      fprintf(srf, "f\n");
   else
      fprintf(srf, "c\n");

   if (numsalts > 0) {
      for (s = 0; s < numsalts-1; s++)
         fprintf(srf, "%s ", saltfilenames[s]);
      fprintf(srf, "%s", saltfilenames[numsalts-1]);
   }
   fprintf(srf, "\n");
   if (numdielectrics > 0) {
      for (s = 0; s < numdielectrics-1; s++)
         fprintf(srf, "%s ", dielectricfilenames[s]);
      fprintf(srf, "%s", dielectricfilenames[numdielectrics-1]);
   }
   fprintf(srf, "\n");
   if (numdielectriccavities > 0) {
      for (s = 0; s < numdielectriccavities-1; s++)
         fprintf(srf, "%s ", dielectriccavityfilenames[s]);
      fprintf(srf, "%s", dielectriccavityfilenames[numdielectriccavities-1]);
   }
   fprintf(srf, "\n");
   if (numsaltcavities > 0) {
      for (s = 0; s < numsaltcavities-1; s++)
         fprintf(srf, "%s ", saltcavityfilenames[s]);
      fprintf(srf, "%s", saltcavityfilenames[numsaltcavities-1]);
   }
   fprintf(srf, "\n");

   for (d = 0; d < numdielectrics; d++) {
      for (s = 0; s < numsalts; s++)
        if (isInside(dielectricpanels[d], numdielectricpanels[d],
                     saltpanels[s], numsaltpanels[s])) {
           fprintf(srf, "%u ", s);
           break;
        }
      if (s == numsalts) {
         if (numsalts > 0)
            printf("Dielectric %u is not inside anything!\n", d);
         fprintf(srf, "999 ");
      }
   }
   fprintf(srf, "\n");

   for (dc = 0; dc < numdielectriccavities; dc++) {
      for (d = 0; d < numdielectrics; d++)
         if (dielectriccode == DIELECTRIC_MSMS_FLAT) {
            if (isInsideExhaustive(dielectriccavitypanels[dc], numdielectriccavitypanels[dc],
                                   dielectricpanels[d], numdielectricpanels[d])) {
               fprintf(srf, "%u ", d);
               break;
            }
         }
         else {
            if (isInside(dielectriccavitypanels[dc], numdielectriccavitypanels[dc],
                         dielectricpanels[d], numdielectricpanels[d])) {
               fprintf(srf, "%u ", d);
               break;
            }
         }
      if (d == numdielectrics) {
         printf("Dielectric cavity %u is not inside anything!\n", dc);
         fprintf(srf, "999 ");
      }
   }
   fprintf(srf, "\n");

   for (sc = 0; sc < numsaltcavities; sc++) {
      for (dc = 0; dc < numdielectriccavities; dc++)
        if (isInside(saltcavitypanels[sc], numsaltcavitypanels[sc],
                     dielectriccavitypanels[dc], numdielectriccavitypanels[dc])) {
           fprintf(srf, "%u ", dc);
           break;
        }
      if (dc == numdielectriccavities) {
         printf("Salt cavity %u is not inside anything!\n", sc);
         fprintf(srf, "999 ");
      }
   }
   fprintf(srf, "\n");

   fclose(srf);
}

int main(int argc, char* argv[]) {
   if (argc != 11) {
      printf("Usage: %s [molecule.pdb/crd/xyzr] [radii.siz] [output.srf] [proberadius] [ionexclusionradius] [dielectricdensity] [saltdensity] [dielectriccode] [saltcode] [prefix]\n", argv[0]);
      return -1;
   }

    if (!strncmp(argv[1] + strlen(argv[1]) - 3, "crd", 3))
       readCRD(argv[1], &numPDBentries, &PDBentries);
    else if (!strncmp(argv[1] + strlen(argv[1]) - 4, "xyzr", 4))
       readXYZR(argv[1], &numPDBentries, &PDBentries);
    else
       readPDB(argv[1], &numPDBentries, &PDBentries);

    if (strncmp(argv[1] + strlen(argv[1]) - 4, "xyzr", 4)) {
       readSIZ(argv[2], &numSIZentries, &SIZentries);
       assignRadiiCharges(PDBentries, numPDBentries, SIZentries, numSIZentries, NULL, 0);
    }

    proberadius = atof(argv[4]);
    ionexclusionradius = atof(argv[5]);
    dielectricdensity = atof(argv[6]);
    saltdensity = atof(argv[7]);
    dielectriccode = atoi(argv[8]);
    saltcode = atoi(argv[9]);

    strcpy(prefix, argv[10]);

    switch (saltcode) {
       case SALT_NONE: break;
       case SALT_MSMS_FLAT: generateSurface_salt_MSMS_flat(); break;
       case SALT_MSMS_CURVED: generateSurface_salt_MSMS_curved(); break;
       case SALT_NETGEN_FLAT: generateSurface_salt_NETGEN_flat(); break;
       case SALT_NETGEN_CURVED: generateSurface_salt_NETGEN_curved(); break;
       default: printf("Unknown salt code %u\n", saltcode); exit(-2); break;
    }

    switch (dielectriccode) {
       case DIELECTRIC_NONE: break;
       case DIELECTRIC_MSMS_FLAT: generateSurface_dielectric_MSMS_flat(); break;
       case DIELECTRIC_MSMS_CURVED: generateSurface_dielectric_MSMS_curved(); break;
       case DIELECTRIC_NETGEN_FLAT: generateSurface_dielectric_NETGEN_flat(); break;
       case DIELECTRIC_NETGEN_CURVED: generateSurface_dielectric_NETGEN_curved(); break;
       default: printf("Unknown dielectric code %u\n", dielectriccode); exit(-2); break;
    }

    unsigned int s, d;

    printf("Surface filenames:\n");

    printf("Salt:");
    for (s = 0; s < numsalts; s++)
       printf(" %s", saltfilenames[s]);
    printf("\n");

    printf("Dielectric:");
    for (s = 0; s < numdielectrics; s++)
       printf(" %s", dielectricfilenames[s]);
    printf("\n");

    printf("Dielectric Cavity:");
    for (s = 0; s < numdielectriccavities; s++)
       printf(" %s", dielectriccavityfilenames[s]);
    printf("\n");

    printf("Salt Cavity:");
    for (s = 0; s < numsaltcavities; s++)
       printf(" %s", saltcavityfilenames[s]);
    printf("\n");

    writeSRF(argv[3]);

    return 0;
}
