#include "meshmaker.h"
#include <stdarg.h>

/* FFTSVDpbeAPI.c */
real ionexclusionradius = 0.0f;
real proberadius = 0.0f;

void error(char* message, ...) {
   va_list ap;

   va_start(ap, message);

   printf("ERROR: ");
   vprintf(message, ap);
   printf("\n");

   va_end(ap);

   exit(-1);
}

void warning(char* message, ...) {
   va_list ap;

   va_start(ap, message);

   printf("WARNING: ");
   vprintf(message, ap);
   printf("\n");

   va_end(ap);
}

void removeWhitespace(char* s) {
   int length = strlen(s);
   char* temp = (char*)calloc(length+1, sizeof(char));
   int scount = 0, tempcount = 0;

   memset(temp, 0, length+1);

   for (scount = 0; scount < length; scount++)
      if (!isspace(s[scount])) {
         temp[tempcount] = s[scount];
         tempcount++;
      }

   strcpy(s, temp);

   free(temp);
}

void readPDB(const char* filename, unsigned int* numPDBentries, PDBentry** PDBentries) {
   FILE* pdbfile = NULL;
   char line[81], number[9];
   int atomcount = 0;

   pdbfile = fopen(filename, "r");

   if (!pdbfile)
      error("Could not open PDB file %s: %s", filename, strerror(errno));

   /* First, read through the file counting the number of ATOM or HETATM lines */

   while (fgets(line, 81, pdbfile)) {
      if (!strncasecmp(line, "ATOM", 4) || !strncasecmp(line, "HETATM", 6))
         atomcount++;
   }

   /* Allocate the memory */

   *PDBentries = (PDBentry*)calloc(atomcount, sizeof(PDBentry));

   if (!(*PDBentries))
      error("Error in memory allocation: %s", strerror(errno));

   /* Rewind the PDB file */

   rewind(pdbfile);

   /* Now read through again and store the entries */

   *numPDBentries = atomcount;

   atomcount = 0;

   while (fgets(line, 81, pdbfile)) {
      if (!strncasecmp(line, "ATOM", 4) || !strncasecmp(line, "HETATM", 6)) {
         strncpy((*PDBentries)[atomcount].record, line, 6);
         strncpy(number, line + 6, 5);
         number[5] = '\0';
         (*PDBentries)[atomcount].atomnumber = atoi(number);
         strncpy((*PDBentries)[atomcount].atomname, line + 12, 4);
         (*PDBentries)[atomcount].alternatelocation = line[16];
         strncpy((*PDBentries)[atomcount].residuename, line + 17, 3);
         (*PDBentries)[atomcount].chain = line[21];
         strncpy(number, line + 22, 4);
         number[4] = '\0';
         (*PDBentries)[atomcount].residuenumber = atoi(number);
         (*PDBentries)[atomcount].residueinsertion = line[26];
         strncpy(number, line + 30, 8);
         number[8] = '\0';
         (*PDBentries)[atomcount].x = atof(number);
         strncpy(number, line + 38, 8);
         number[8] = '\0';
         (*PDBentries)[atomcount].y = atof(number);
         strncpy(number, line + 46, 8);
         number[8] = '\0';
         (*PDBentries)[atomcount].z = atof(number);
         strncpy(number, line + 54, 6);
         number[6] = '\0';
         (*PDBentries)[atomcount].occupancy = atof(number);
         strncpy(number, line + 60, 6);
         number[6] = '\0';
         (*PDBentries)[atomcount].temperature = atof(number);
         strncpy(number, line + 67, 3);
         number[3] = '\0';
         (*PDBentries)[atomcount].footnotenumber = atof(number);

         removeWhitespace((*PDBentries)[atomcount].record);
         removeWhitespace((*PDBentries)[atomcount].atomname);
         removeWhitespace((*PDBentries)[atomcount].residuename);

         atomcount++;
      }
   }

   fclose(pdbfile);
}

void readCRD(const char* filename, unsigned int* numPDBentries, PDBentry** PDBentries) {
   FILE* crdfile = NULL;
   char line[81], number[11];
   int atomcount = 0;

   crdfile = fopen(filename, "r");

   if (!crdfile)
      error("Could not open CRD file %s: %s", filename, strerror(errno));

   /* First, read through the header lines */

   while (fgets(line, 81, crdfile))
      if (line[0] != '*')
         break;

   atomcount = atoi(line);

   /* Allocate the memory */

   *PDBentries = (PDBentry*)calloc(atomcount, sizeof(PDBentry));

   if (!(*PDBentries))
      error("Error in memory allocation: %s", strerror(errno));

   /* Now read through again and store the entries */

   *numPDBentries = atomcount;

   atomcount = 0;

   while (fgets(line, 81, crdfile)) {
      strncpy(number, line, 5);
      number[5] = '\0';
      (*PDBentries)[atomcount].atomnumber = atoi(number);
      strncpy(number, line + 5, 5);
      number[5] = '\0';
      (*PDBentries)[atomcount].residuenumber = atoi(number);
      strncpy((*PDBentries)[atomcount].residuename, line + 11, 3);
      strncpy((*PDBentries)[atomcount].atomname, line + 16, 4);
      strncpy(number, line + 20, 10);
      number[10] = '\0';
      (*PDBentries)[atomcount].x = atof(number);
      strncpy(number, line + 30, 10);
      number[10] = '\0';
      (*PDBentries)[atomcount].y = atof(number);
      strncpy(number, line + 40, 10);
      number[10] = '\0';
      (*PDBentries)[atomcount].z = atof(number);
      (*PDBentries)[atomcount].chain = line[51];
      strncpy(number, line + 61, 10);
      number[10] = '\0';
      (*PDBentries)[atomcount].temperature = atof(number);

      removeWhitespace((*PDBentries)[atomcount].atomname);
      removeWhitespace((*PDBentries)[atomcount].residuename);

      atomcount++;
   }

   fclose(crdfile);
}

void readXYZR(const char* filename, unsigned int* numPDBentries, PDBentry** PDBentries) {
   FILE* xyzrfile = NULL;
   char line[81], number[11];
   int atomcount = 0;

   xyzrfile = fopen(filename, "r");

   if (!xyzrfile)
      error("Could not open XYZR file %s: %s", filename, strerror(errno));

   /* First, read through the file counting the number of lines */

   while (fgets(line, 81, xyzrfile))
      atomcount++;

   /* Allocate the memory */

   *PDBentries = (PDBentry*)calloc(atomcount, sizeof(PDBentry));

   if (!(*PDBentries))
      error("Error in memory allocation: %s", strerror(errno));

   /* Rewind the PDB file */

   rewind(xyzrfile);

   /* Now read through again and store the entries */

   *numPDBentries = atomcount;

   atomcount = 0;

   while (fgets(line, 81, xyzrfile)) {
      sscanf(line, "%lf %lf %lf %lf", &(*PDBentries)[atomcount].x, &(*PDBentries)[atomcount].y, &(*PDBentries)[atomcount].z, &(*PDBentries)[atomcount].radius);
      atomcount++;
   }

   fclose(xyzrfile);
}

void readSIZ(const char* filename, unsigned int* numSIZentries, SIZentry** SIZentries) {
   FILE* sizfile = NULL;
   char line[81], number[9];
   int radiicount = 0;

   sizfile = fopen(filename, "r");

   if (!sizfile)
      error("Could not open SIZ file %s: %s", filename, strerror(errno));

   /* First, read through the file counting the number of non comment lines */

   fgets(line, 81, sizfile);  /* first line is junk */

   while (fgets(line, 81, sizfile)) {
      if ((line[0] != '!') && (strlen(line) > 0))
         radiicount++;
   }

   /* Allocate the memory */

   *SIZentries = (SIZentry*)calloc(radiicount, sizeof(SIZentry));

   if (!(*SIZentries))
      error("Error in memory allocation: %s", strerror(errno));

   /* Rewind the SIZ file */

   rewind(sizfile);

   /* Now read through again and store the entries */

   *numSIZentries = radiicount;

   radiicount = 0;

   fgets(line, 81, sizfile);  /* first line is junk */

   while (fgets(line, 81, sizfile)) {
      if ((line[0] != '!') && (strlen(line) > 0)) {
         strncpy((*SIZentries)[radiicount].atomlabel, line, 6);
         strncpy((*SIZentries)[radiicount].residuelabel, line + 6, 3);
         strncpy(number, line + 9, 8);
         number[8] = '\0';
         (*SIZentries)[radiicount].radius = atof(number);

         removeWhitespace((*SIZentries)[radiicount].atomlabel);
         removeWhitespace((*SIZentries)[radiicount].residuelabel);

         radiicount++;
      }
   }

   fclose(sizfile);
}

void readCRG(const char* filename, unsigned int* numCRGentries, CRGentry** CRGentries) {
#ifndef SCATTER
   FILE* crgfile = NULL;
   char line[81], number[9];
   int chargecount = 0;

   crgfile = fopen(filename, "r");

   if (!crgfile)
      error("Could not open CRG file %s: %s", filename, strerror(errno));

   /* First, read through the file counting the number of non comment lines */

   fgets(line, 81, crgfile);  /* first line is junk */

   while (fgets(line, 81, crgfile)) {
      if ((line[0] != '!') && (strlen(line) > 0))
         chargecount++;
   }

   /* Allocate the memory */
	if (!(*CRGentries))
	  free(*CRGentries);
	
   *CRGentries = (CRGentry*)calloc(chargecount, sizeof(CRGentry));

   if (!(*CRGentries))
      error("Error in memory allocation: %s", strerror(errno));

   /* Rewind the CRG file */

   rewind(crgfile);

   /* Now read through again and store the entries */

   *numCRGentries = chargecount;

   chargecount = 0;

   fgets(line, 81, crgfile);  /* first line is junk */

   while (fgets(line, 81, crgfile)) {
      if ((line[0] != '!') && (strlen(line) > 0)) {
         strncpy((*CRGentries)[chargecount].atomlabel, line, 6);
         strncpy((*CRGentries)[chargecount].residuelabel, line + 6, 3);
         strncpy(number, line + 9, 4);
         number[4] = '\0';
         (*CRGentries)[chargecount].residuenumber = atoi(number);
         (*CRGentries)[chargecount].chain = line[13];
         strncpy(number, line + 14, 8);
         number[8] = '\0';
         (*CRGentries)[chargecount].charge = atof(number);

         removeWhitespace((*CRGentries)[chargecount].atomlabel);
         removeWhitespace((*CRGentries)[chargecount].residuelabel);

         chargecount++;
      }
   }

   fclose(crgfile);
#else
	printf("Skipping readCRG because scattering calculation enabled!\n");
	*numCRGentries = 0;
	*CRGentries = NULL;
#endif
}

void writeCRG(char* filename, PDBentry *PDBentries, unsigned int numPDBentries)
{
  FILE* crgfile = NULL;
  char line[25];
  crgfile = fopen(filename, "w");
  
  //line = malloc(sizeof(char)*100);
  
  if(!crgfile)
    error("Could not open CRG file %s for writing: %s", filename, strerror(errno));

  fprintf(crgfile, "aaaaaarrrnnnncqqqqqqqq\n");
  



  for(int i=0; i<numPDBentries; ++i) {
    for(int i=0; i<25; ++i)
      line[i] = ' ';
    line[24] = '\0';
    line[23] = '\n';
    strncpy(line, PDBentries[i].atomname, 6);
    int size = strlen(PDBentries[i].atomname);
    for(int j=size; j<6; ++j)
      line[j] = ' ';
    strncpy(line+6, PDBentries[i].residuename, 3);
    size = strlen(PDBentries[i].residuename);
    printf("RESSIZE: %d\n", size);
    for(int j=size+6; j<9; ++j)
      line[j] = ' ';
    snprintf(line+14, sizeof(char)*9, "% 7.5f", PDBentries[i].charge);
    fputs(line, crgfile);
    fputs("\n", crgfile);
  }
  
  fclose(crgfile);

}

void assignRadiiCharges(PDBentry* PDBentries, unsigned int numPDBentries, SIZentry* SIZentries, unsigned int numSIZentries, CRGentry* CRGentries, unsigned int numCRGentries) {
#ifndef SCATTER
   unsigned int a, r, c, matchlevel;

   for (a = 0; a < numPDBentries; a++) {
      matchlevel = 0;

      for (r = 0; r < numSIZentries; r++) {
         if ((!strcmp(PDBentries[a].atomname, SIZentries[r].atomlabel)) && (!strcmp(PDBentries[a].residuename, SIZentries[r].residuelabel))) {
            PDBentries[a].radius = SIZentries[r].radius;
            matchlevel = 3;
            break;  /* A complete match, so we can stop here */
         }
         else if ((strlen(SIZentries[r].residuelabel) == 0) && (!strncmp(PDBentries[a].atomname, SIZentries[r].atomlabel, strlen(SIZentries[r].atomlabel)))) {
            if (matchlevel < 2) {
               PDBentries[a].radius = SIZentries[r].radius;
               matchlevel = 2;
            }
         }
         else if ((strlen(SIZentries[r].residuelabel) == 0) && (PDBentries[a].atomname[0] == SIZentries[r].atomlabel[0]) && (strlen(SIZentries[r].atomlabel) == 1)) {
            if (matchlevel < 1) {
               PDBentries[a].radius = SIZentries[r].radius;
               matchlevel = 1;
            }
         }
      }

      if (matchlevel == 0) {
         warning("Atom %d %s %s %c has no radius entry, defaulting to 0.0",
                 PDBentries[a].atomnumber, PDBentries[a].atomname,
                 PDBentries[a].residuename, PDBentries[a].chain);

         PDBentries[a].radius = 0.0f;
      }
      else if (matchlevel == 1)
         warning("Atom %d %s %s %c using default element radius of %f",
                 PDBentries[a].atomnumber, PDBentries[a].atomname,
                 PDBentries[a].residuename, PDBentries[a].chain,
                 PDBentries[a].radius);

      matchlevel = 0;

      for (c = 0; c < numCRGentries; c++) {
	if ((CRGentries[c].chain == PDBentries[a].chain) && (CRGentries[c].residuenumber == PDBentries[a].residuenumber) && (!strcmp(PDBentries[a].residuename, CRGentries[c].residuelabel)) && (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel))) {
            PDBentries[a].charge = CRGentries[c].charge;
            matchlevel = 8;
            break;  /* A complete match, so we can stop here */
         }
         else if ((CRGentries[c].residuenumber == PDBentries[a].residuenumber) && (!strcmp(PDBentries[a].residuename, CRGentries[c].residuelabel)) && (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel))) {
            if (matchlevel < 7) {
               PDBentries[a].charge = CRGentries[c].charge;
               matchlevel = 7;
            }
         }
         else if ((CRGentries[c].chain == PDBentries[a].chain) && (!strcmp(PDBentries[a].residuename, CRGentries[c].residuelabel)) && (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel))) {
            if (matchlevel < 6) {
               PDBentries[a].charge = CRGentries[c].charge;
               matchlevel = 6;
            }
         }
         else if ((!strcmp(PDBentries[a].residuename, CRGentries[c].residuelabel)) && (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel))) {
            if (matchlevel < 5) {
               PDBentries[a].charge = CRGentries[c].charge;
               matchlevel = 5;
            }
         }
         else if ((CRGentries[c].chain == PDBentries[a].chain) && (CRGentries[c].residuenumber == PDBentries[a].residuenumber) && (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel))) {
            if (matchlevel < 4) {
               PDBentries[a].charge = CRGentries[c].charge;
               matchlevel = 4;
            }
         }
         else if ((CRGentries[c].residuenumber == PDBentries[a].residuenumber) && (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel))) {
            if (matchlevel < 3) {
               PDBentries[a].charge = CRGentries[c].charge;
               matchlevel = 3;
            }
         }
         else if ((CRGentries[c].chain == PDBentries[a].chain) && (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel))) {
            if (matchlevel < 2) {
               PDBentries[a].charge = CRGentries[c].charge;
               matchlevel = 2;
            }
         }
         else if (!strcmp(PDBentries[a].atomname, CRGentries[c].atomlabel)) {
            if (matchlevel < 1) {
               PDBentries[a].charge = CRGentries[c].charge;
               matchlevel = 1;
            }
         }
      }

      if ((matchlevel == 0) && (numCRGentries > 0)) {
         warning("Atom %d %s %s %c has no charge entry, defaulting to 0.0",
                 PDBentries[a].atomnumber, PDBentries[a].atomname,
                 PDBentries[a].residuename, PDBentries[a].chain);

         PDBentries[a].charge = 0.0f;
      }
   }
#else
	printf("Skipping assignRadiiCharges because scattering enabled.\n");
#endif 
}
