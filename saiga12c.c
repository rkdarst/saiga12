/* Richard Darst, 2008
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#define CHECK printf(".\n");

#include "SFMT.h"

#define S12_EMPTYSITE (-1)
#define S12_TYPE_ANY (-2)
#define S12_ENERGY_BM (1)
#define S12_ENERGY_ZERO (2)
#define S12_ENERGY_CTCC (3)
#define S12_ENERGY_SPM (4)
#define S12_ENERGY_TPM (5)
#define S12_ENERGY_BMnotzero (10)
#define S12_ENERGY_BMimmobile1 (11)
#define S12_ENERGY_BMimmobile1b (12)
#define S12_CYCLE_MC (1)
#define S12_CYCLE_KA (2)
#define S12_CYCLE_FA (3)
#define S12_CYCLE_CTCC (4)
#define S12_CYCLE_EAST (5)
#define S12_CYCLE_SPIRAL (6)
#define S12_CYCLE_SPINMC (7)
#define S12_CYCLE_CTCCclassic (10)
#define S12_FLAG_VIB_ENABLED (1)
#define S12_FLAG_DOSIN (2)
#define S12_FLAG_FROZEN (4)
#define S12_FLAG_SELECTED (8)
#define S12_FLAG_KA_SOFT (16)

#define S12_FLAG_INCOMPAT (S12_FLAG_FROZEN)


//These all used to be globally defined ints as flags.
#ifndef debug
#define debug 0
#endif
// Check for consistency errors and print if found
#ifndef errorcheck
#define errorcheck 0
#endif
// Log event-driven dynamics code
#ifndef debugedd
#define debugedd 0
#endif

struct SimData {
  double beta;
  int N;
  int *ntype;
  int ntypeMax;
  //int NMax;  // lattSize is NMax
  double hardness;
  int cycleMode;   // see python for definition
  int energyMode;  // see python for definition
  int flags;
  int (*callback)(struct SimData *SD);
  void *S;

  double uVTchempotential;
  int inserttype;  // must be generalized later.
  int widominserttype;  // must be generalized later.
  int inserttypes_n;    // not used
  double *inserttypes_prob; // when multiple insert types, this is cum prob.
  int *inserttypes_type; // and this is corresponding type.
  double *inserttypes_plookup; // lookup from type->prob
  double *inserttypes_mulookup; // lookup from type->mu

  int    *MLL;    // Move Lookup List, for event-driven dynamics
  int    *MLLr;   // Move Lookup List, reverse
  int    MLLlen;
  int    *MLL_down;    // Move Lookup List, for FA down spins.
  int    MLLlen_down;
  double MLLextraTime;

  int lattSize;
  int *lattsite;
  int *conn;
  int *connN;
  int connMax;
  int *nneighbors;
  int *atomtype;
  int *atompos;
  int *persist;
  int *orient;
  int *frozen;
  int *selected;

  double cumProbAdd;
  double cumProbDel;
} ;

struct LList {
  int n;
  int elem[30];
} ;
#define LListN 30

/* genrand_real1  [0,1]
 * genrand_real2  [0,1)
 * genrand_real3  (0,0)
 */


void ctest(struct SimData *SD);

/* Given lattice point `pos`, how many adjacent atoms does it have?
 */
inline int neighbors_pos(struct SimData *SD, int pos) {
  int neighbors = 0;
  int i;
  for (i=0 ; i<SD->connN[pos] ; i++) {
    int neighpos = SD->conn[SD->connMax*pos + i];
    if (SD->lattsite[neighpos] != S12_EMPTYSITE) {
      neighbors += 1;
    }
  }
  return(neighbors);
}

/*  Use the getInsertTypes returns a random type to try inserting in
    multicomponent grand canonical simulations.
 */
//#define getInsertType(x) (x)->inserttype;
inline int getInsertType(struct SimData *SD) {
  /* When running a mixture, this selects of several different
   * particle types to insert.
   * - We have arrays inserttypes_prob and inserttypes_type, which hold:
   *   a) the upper bound cumulative prob function
   *   b) the type corresponding to that upper bound CPDF.
   */
  if (SD->inserttype != S12_EMPTYSITE)
    return(SD->inserttype);
  else {
    double ran = genrand_real2();  /*[0,1)*/
    int i=0;
    while (1) {
      if (SD->inserttypes_prob[i] >= ran)
	return(SD->inserttypes_type[i]);
      i++;
    }
  }
}

/*  C function to print out the lattice.  Not really used except for
 *  debugging from C-- use the python ones instead.
 */
void printLatticeC(struct SimData *SD) {
  int pos;
  for (pos=0 ; pos<SD->lattSize ; pos++)
    printf("%d ", SD->lattsite[pos]);
  printf("\n");
}

/*  Little "macro" to return atomtype at a position.
 */
inline int atomType(struct SimData *SD, int pos) {
  //return(SD->lattsite[pos]);
  return(SD->atomtype[SD->lattsite[pos]]);
}
/* Little macro for finding the reverse connection.
 */
inline int reverseConnection(struct SimData *SD, int i_conn) {
  int connMax = SD->connMax;
  return ( (i_conn+(connMax/2)) % connMax );
}


inline void addParticle(struct SimData *SD, int pos, int type) {
  if (errorcheck) if (SD->lattsite[pos] != S12_EMPTYSITE) { 
      printf("error: inserting atom at not empty site location: %d\n", pos);
      exit(61); }
  if (errorcheck) if (type > SD->ntypeMax || type < 0) {
    printf("Inserting particle of type>ntypeMax or type<zero: %d\n", type);
    exit(60);
  }

  ////SD->lattsite[pos] = type;
  SD->lattsite[pos] = SD->N;     // atomnumberings start at zero.
  int i;
  for (i=0 ; i<SD->connN[pos] ; i++) {
    int neighpos = SD->conn[SD->connMax*pos + i];
    SD->nneighbors[neighpos] ++;
  }
  SD->atomtype[SD->N] = type;
  SD->atompos[SD->N] = pos;
  SD->ntype[type]++;
  if (SD->orient != NULL) {
    int conni = SD->connN[pos] * genrand_real2();
    SD->orient[pos] = conni;
  }
  SD->N++;
}
inline void delParticle(struct SimData *SD, int pos) {
  if (errorcheck) if (SD->lattsite[pos] == S12_EMPTYSITE) {
      printf("error: removing atom at empty site location: %d\n", pos);
      exit(62); }
  int i;
  for (i=0 ; i<SD->connN[pos] ; i++) {
    int neighpos = SD->conn[SD->connMax*pos + i];
    SD->nneighbors[neighpos] --;
  }
  SD->ntype[atomType(SD, pos)]--;
  // simply remove it if it's the last one
  int atomnum = SD->lattsite[pos];
  if (atomnum == SD->N-1) {  // easy case: just remove
    //atomType
    SD->atomtype[atomnum] = S12_EMPTYSITE;  // reset it
    SD->atompos[atomnum] = S12_EMPTYSITE;  // reset it
  }
  else {
    // move the N-1 atom to atomnum.
    SD->atomtype[atomnum] = SD->atomtype[SD->N-1];
    SD->atompos[atomnum] = SD->atompos[SD->N-1];
    SD->lattsite[SD->atompos[SD->N-1]] = atomnum;

    SD->atomtype[SD->N-1] = S12_EMPTYSITE;  // reset it
    SD->atompos[SD->N-1] = S12_EMPTYSITE;  // reset it
  }
  SD->lattsite[pos] = S12_EMPTYSITE;
  if (SD->orient != NULL) {
    SD->orient[pos] = S12_EMPTYSITE ;
  }
  SD->N--;
}
inline void moveParticle(struct SimData *SD, int oldpos, int newpos) {
  if (errorcheck) if (SD->lattsite[newpos] != S12_EMPTYSITE) {
      printf("error: move atom to not empty site location: %d\n", newpos);
      exit(63); }
  if (errorcheck) if (SD->lattsite[oldpos] == S12_EMPTYSITE) {
      printf("error: move atom from empty site location: %d\n", oldpos);
      exit(64); }
  SD->lattsite[newpos] = SD->lattsite[oldpos];
  SD->lattsite[oldpos] = S12_EMPTYSITE;
  SD->atompos[SD->lattsite[newpos]] = newpos;
  int i;
  for (i=0 ; i<SD->connN[oldpos] ; i++) {
    int neighpos = SD->conn[SD->connMax*oldpos + i];
    SD->nneighbors[neighpos] --;
  }
  for (i=0 ; i<SD->connN[newpos] ; i++) {
    int neighpos = SD->conn[SD->connMax*newpos + i];
    SD->nneighbors[neighpos] ++;
  }
  if (SD->orient != NULL) {
    if(errorcheck) {
      if(SD->conn[SD->connMax*oldpos+SD->orient[oldpos]] != newpos) {
	printf("moved in direction we aren't oriented in."); exit (79);
      }
    }
    // what is our new orientation?
    int oldOrient = SD->orient[oldpos];
    int newOrient = (oldOrient+(SD->connMax/2)) % SD->connMax;
    // which direction should it now be facing?
    SD->orient[newpos] = newOrient;
    SD->orient[oldpos] = S12_EMPTYSITE ;
  }
}
inline void moveParticleShort(struct SimData *SD, int oldpos, int newpos) {
  /* Doesn't errorcheck - used in CTCCclassic dynamics.
   */
  SD->lattsite[newpos] = SD->lattsite[oldpos];
  SD->lattsite[oldpos] = S12_EMPTYSITE;
  SD->atompos[SD->lattsite[newpos]] = newpos;
  int i;
  for (i=0 ; i<SD->connN[oldpos] ; i++) {
    int neighpos = SD->conn[SD->connMax*oldpos + i];
    SD->nneighbors[neighpos] --;
  }
  for (i=0 ; i<SD->connN[newpos] ; i++) {
    int neighpos = SD->conn[SD->connMax*newpos + i];
    SD->nneighbors[neighpos] ++;
  }
}
void randomizeSystem(struct SimData *SD, int type, double density) {
  int pos, n;
  for (n=SD->N-1 ; n >= 0 ; n--) {
    delParticle(SD, SD->atompos[n]);
  }
  int numParticles = SD->lattSize * density;
  for (n=0 ; n<numParticles ; n++) {
    do {
      pos = SD->lattSize * genrand_real2();
    } while ( SD->lattsite[pos] != S12_EMPTYSITE );
    addParticle(SD, pos, type);
  }
}


void loadStateFromSave(struct SimData *SD) {
  /* This function re-initilized the various arrays, from just the
   * lattsite and type arrays.
   */
  int n, pos;
  for (n=0 ; n<SD->N ; n++) {
    SD->lattsite[SD->atompos[n]] = n;
    SD->ntype[SD->atomtype[n]] += 1;
  }
  for (pos=0 ; pos<SD->lattSize ; pos++) {
    SD->nneighbors[pos] = neighbors_pos(SD, pos) ;
  }
}

/* How to use these energy functions ?
 * energy_pos is the main energy function.
 * energy_posLocal is the energy due to one lattice site.  It is only
 *           important when there are "induced 3-body interactions"
 * energy is the sum of all energies in the system.
 * 
 * How to add a new energy function:
 * a) add a energy_XXX.c file in ccode/ directory.
 *    add an #include "ccode/energy_XXX.c" below to include it.
 * b) add S12_ENERGY_XXX constant in the top of saiga12c.c AND common.py
 * c) add the energy setting in setEnergyMode in main.py.  Follow the
 *    examples.
 * d) Most likely, your energy function is two-body.  If so, do this:
 *    d1) you need only one energy function, call it 
 *        energyXXX_pos(sturct SimData *SD, int pos)
 *        in the ccode/ file.  Make it an inline function
 *    d2) add the hooks in energy_posLocal and energy_pos to call the
 *        function above.  You don't need to make energy_pos and
 *        energy_posLocal call different things since your function is
 *        only two-body.
 * e) saiga12 extensively uses inline functions in critical portions.
 *    It would probably go faster if un-inlined the energy functions you
 *    were not using -- so remove the inlines from ccode/energy_bm.c.
 */
#include "ccode/energy_bm.c"
#include "ccode/energy_ctcc.c"
#include "ccode/energy_squareplaquette.c"
#include "ccode/energy_triangularplaquette.c"

inline double energy_posLocal(struct SimData *SD, int pos) {
  /*  Energy of a particle at one particular position.  This function
   *  should only be used instead of energy_pos when there are three
   *  body interactions, and you want to isolate which body is causing
   *  a problem.  In short, don't use it unless you know what you are
   *  doing.
   */
  if (SD->energyMode == S12_ENERGY_BM)
    return energyBM_posLocal(SD, pos);
  else if (SD->energyMode == S12_ENERGY_ZERO)
    return(0.);
  else if (SD->energyMode == S12_ENERGY_CTCC)
    return energyCTCC_posLocal(SD, pos);
  else if (SD->energyMode == S12_ENERGY_BMnotzero)
    return energyBMnotzero_posLocal(SD, pos);
  else if (SD->energyMode == S12_ENERGY_BMimmobile1)
    return energyBMimmobile1_posLocal(SD, pos);
  else if (SD->energyMode == S12_ENERGY_BMimmobile1b)
    return energyBMimmobile1b_posLocal(SD, pos);
  else if (SD->energyMode == S12_ENERGY_SPM)
    return energySPM_posLocal(SD, pos);
  else if (SD->energyMode == S12_ENERGY_TPM)
    return energyTPM_posLocal(SD, pos);
  // add new energy modes above this line.
  else {
    if (errorcheck) {
      printf("invalid evergy mode: %d", SD->energyMode);
      exit(5); }
    else
      return(-1.);
  }
}

inline double energy_pos(struct SimData *SD, int pos) {
  /* This should be your primary energy function for a certain lattice site.
   */
  if (SD->energyMode == S12_ENERGY_BM)
    return energyBM_pos(SD, pos);
  else if (SD->energyMode == S12_ENERGY_ZERO)
    return(0.);
  else if (SD->energyMode == S12_ENERGY_CTCC)
    return energyCTCC_pos(SD, pos);
  else if (SD->energyMode == S12_ENERGY_BMnotzero)
    return energyBMnotzero_pos(SD, pos);
  else if (SD->energyMode == S12_ENERGY_BMimmobile1)
    return energyBMimmobile1_pos(SD, pos);
  else if (SD->energyMode == S12_ENERGY_BMimmobile1b)
    return energyBMimmobile1b_pos(SD, pos);
  else if (SD->energyMode == S12_ENERGY_SPM)
    return energySPM_pos(SD, pos);
  else if (SD->energyMode == S12_ENERGY_TPM)
    return energyTPM_pos(SD, pos);
  // add new energy modes above this line.
  else {
    if (errorcheck) {
      printf("invalid evergy mode: %d", SD->energyMode);
      exit(6); }
    else
      return(-1.);
  }
}

double energy(struct SimData *SD) {
  /* Total energy of the system.  Note that for two-body symmetric
   * potentials, this does not divide by two to eliminate
   * double-counting.
   */
  int pos;
  double energy=0;
  for (pos=0 ; pos<SD->lattSize ; pos++) {
    energy += energy_posLocal(SD, pos);
  }
  return energy;
}



/* These functions make a simple list of numbers.  They are used to
 * optimize later code, so that when a site is processed, you add it
 * to the list with LlistAdd.  Then, before processing again, you
 * check with LlistLooup to see if it has already been processed.
 *
 * So basically it as a data structure where you can a) add things b)
 * see if an element is already in it c) zero it (letting
 * llist.len=0).
 */
inline void LlistAdd(struct LList *llist, int elem) {
  if (llist->n == LListN) {
    return ; // We are already full, silently ignore adding it.
  }
  llist->elem[llist->n] = elem ;
  llist->n ++;
}
inline int LlistLookup(struct LList *llist, int elem) {
  int i = llist->n;
  //printf("Doing LlistLookup: n:%d\n", llist->n);
  while (i--) {
    //printf(" ..LlistLookup: %d %d\n", llist->n, i);
    if (llist->elem[i] == elem)
      return 1;
  }
  return 0;
}




inline void addToMLL(struct SimData *SD, int pos, int conni) {
  int lookup = pos * SD->connMax + conni;
  if (errorcheck && (SD->MLLr[lookup] != S12_EMPTYSITE)) {
    printf("error rhoeacurk %d %d \n", pos, conni);
  }
  SD->MLL[SD->MLLlen] = lookup;
  SD->MLLr[lookup] = SD->MLLlen;
  SD->MLLlen += 1;
}
inline void removeFromMLL(struct SimData *SD, int pos, int conni) {
  int moveIndex = pos * SD->connMax + conni;
  if (errorcheck && (SD->MLLr[moveIndex] == S12_EMPTYSITE)) {
    printf("error torknotur %d %d \n", pos, conni);
  }
  int mllLocation = SD->MLLr[moveIndex];
  if (mllLocation == SD->MLLlen-1 ) {
    SD->MLL[mllLocation] = -1;
    SD->MLLr[moveIndex] = -1;
  } else {
    SD->MLL[mllLocation] = SD->MLL[SD->MLLlen-1];
    SD->MLLr[SD->MLL[SD->MLLlen-1]] = mllLocation;

    SD->MLLr[moveIndex] = -1;
    SD->MLL[SD->MLLlen-1] = -1;
  }
  SD->MLLlen -= 1 ;
}

inline int isInMLL(struct SimData *SD, int pos, int conni) {
  if (SD->MLLr[pos*SD->connMax + conni] == S12_EMPTYSITE)
    return 0;
  else
    return 1;
}
inline void ensureInMLL(struct SimData *SD, int pos, int conni) {
  if (SD->MLLr[pos*SD->connMax + conni] == S12_EMPTYSITE)
    addToMLL(SD, pos, conni);
}
inline void ensureNotInMLL(struct SimData *SD, int pos, int conni) {
  if (SD->MLLr[pos*SD->connMax + conni] != S12_EMPTYSITE)
    removeFromMLL(SD, pos, conni);
}
inline void ensureInMLLIf(int isAllowed, struct SimData *SD,int pos,int conni){
  if (isAllowed)
    ensureInMLL(SD, pos, conni);
  else
    ensureNotInMLL(SD, pos, conni);
}
#include "ccode/spinedd.c"


int cycleMC(struct SimData *SD, double n);
int cycleKA(struct SimData *SD, double n);
int cycleFA(struct SimData *SD, double n);
int cycleCTCC(struct SimData *SD, double n);
int cycleCTCCclassic(struct SimData *SD, double n);
inline int cycleKA_translate(struct SimData *SD);
int cycleSpinMC(struct SimData *SD, double n);


int cycle(struct SimData *SD, double n) {
  // If there is no cycle mode set (defaults to empty zero), we should
  // have an error.
  if (SD->cycleMode == S12_CYCLE_MC)
    return cycleMC(SD, n);
  else if (SD->cycleMode == S12_CYCLE_KA)
    return cycleKA(SD, n);
  else if (SD->cycleMode == S12_CYCLE_FA)
    return cycleFA(SD, n);
  else if (SD->cycleMode == S12_CYCLE_CTCC)
    return cycleCTCC(SD, n);
  else if (SD->cycleMode == S12_CYCLE_CTCCclassic)
    return cycleCTCCclassic(SD, n);
  else if (SD->cycleMode == S12_CYCLE_SPINMC)
    return cycleSpinMC(SD, n);
  else {
    printf("Cycle mode not set: %d", SD->cycleMode);
    exit(49);
  }
}

#include "ccode/montecarlo.c"

inline void EddBM_updateLatPos(struct SimData *SD, int pos);
#include "ccode/birolimezard.c"

inline void EddKA_updateLatPos(struct SimData *SD, int pos);
#include "ccode/kobandersen.c"

inline void EddFA_updateLatPos(struct SimData *SD, int pos);
#include "ccode/fredricksonandersen.c"

inline void EddCTCC_updateLatPos(struct SimData *SD, int pos);
#include "ccode/ctcc.c"
#include "ccode/ctccclassic.c"

inline void EddEast_updateLatPos(struct SimData *SD, int pos);
#include "ccode/east.c"

inline void EddSpiral_updateLatPos(struct SimData *SD, int pos);
#include "ccode/spiral.c"

#include "ccode/spinmontecarlo.c"




double chempotential(struct SimData *SD, int inserttype) {
  /*  Use widom insertion to calculate the chemical potential of a
   *  certain type.  This tries inserting a particle at *every*
   *  lattice location, for best averaging.
   */
  // this variable will hold the sum of all exp(-beta deltaU)
  double totalsum=0;
  //int inserttype = SD->widominserttype;
  //printf("insert type: %d\n", inserttype);

  int pos;
  // The for loop takes the average over all possible insertion
  // positions.
  for(pos=0 ; pos<SD->lattSize ; pos++) {
    if (SD->lattsite[pos] != S12_EMPTYSITE) {
      // energy becomes inf, exp() becomes zero, no addition to sum,
      // so continue
      continue;
    }
    double Eold = energy_pos(SD, pos);
    addParticle(SD, pos, inserttype);
    double Enew = energy_pos(SD, pos);
    delParticle(SD, pos);
    //printf("%f %f\n", Eold, Enew);

    totalsum += exp(-SD->beta * (Enew-Eold));
  }
  //printf("%f\n", totalsum);

  //double chempotential = - log(totalsum / SD->lattSize) / SD->beta; // A
  //double chempotential = - log(totalsum / (SD->N+1)) / SD->beta; // B
  double chempotential = 
     - log(totalsum / (SD->ntype[inserttype])) / SD->beta; // B per type
  return(chempotential);
}
double chempotential_innersum(struct SimData *SD, int inserttype) {
  double totalsum=0;
  int pos;
  for(pos=0 ; pos<SD->lattSize ; pos++) {
    if (SD->lattsite[pos] != S12_EMPTYSITE) {
      continue;
    }
    double Eold = energy_pos(SD, pos);
    addParticle(SD, pos, inserttype);
    double Enew = energy_pos(SD, pos);
    delParticle(SD, pos);
    totalsum += exp(-SD->beta * (Enew-Eold));
  }
  return(totalsum / (SD->ntype[inserttype]));
}



double calc_structfact(struct SimData *SD1, struct SimData *SD2,
		       int *atomlist, int N,
		       double *kvecs, int Nk,
		       double *cords, double *cords2,
		       double *lattShape, int nDim,
		       double *Skresult,
		       double *SkArrayByAtom,
		       int flags) {
  int n1, nk;
  double totalsum = 0;
  int pos1, pos2;
  //double r1[3];
  double dr[5];  // max number of dims
  int print=0;

  for (n1=0 ; n1 < N ; n1++) {
    //printf("%p %p\n", SD1, SD2);
    if (print) printf("=== n: %d ===\n", n1);
    int atomn = atomlist[n1];
    pos1 = SD1->atompos[atomn];
    pos2 = SD2->atompos[atomn];
    if ((pos1 == pos2) && (! (flags & S12_FLAG_VIB_ENABLED))) {
      // This simplifying branch is only valid if we do NOT have vibrations.
      totalsum += Nk;
      SkArrayByAtom[n1] += Nk;
      for(nk=0 ; nk < Nk ; nk++)
	Skresult[nk] += 1.;
      continue;
    }

    int d;
    for (d=0 ; d<nDim ; d++) {
      //printf("coord1, coord2:  %f %f \n",
      //     (double)cords[nDim*pos1 + d], (double)cords2[nDim*pos2 + d]);
      //printf(" shape: %f\n", (double)lattShape[d]);
      dr[d] =  cords[nDim*pos1 + d] - cords2[nDim*pos2 + d];
      dr[d] -= (floor(dr[d]/lattShape[d] + .5)) *lattShape[d];
    }
    if(print) printf("  dr: %f %f %f (%d)\n", dr[0], dr[1], dr[2], n1);
    
	
    for(nk=0 ; nk < Nk ; nk++) {
      double dot=0;
      for (d=0 ; d<nDim ; d++) 
	dot += dr[d]*kvecs[nDim*nk+d];
      
      double x = cos(dot);
      if(print) printf("  dot:%f  x:%f   [%f %f %f]\n", dot, x,
		       kvecs[3*nk+0], kvecs[3*nk+1], kvecs[3*nk+2]);
      totalsum += x;
      SkArrayByAtom[n1] += x;
      Skresult[nk] += x;
    }
    //}
  }
  //printf("totalsum: %f ", totalsum);
  return(totalsum);
}



int fourpoint(struct SimData *SD1, struct SimData *SD2,
		 double *cords, double *cords2,
		 double *lattShape, int nDim, int type,
		 double *kvecs, int Nk,
		 double *q,
		 double *A_, double *B_, double *C_,
		 double *Skresult, double *SkArrayByAtom,
		 int flags) {
  int n_part, n_kvec;
  int pos1, pos2;
  double dr[5];  // max number of dims
  int print=0;
  int d;
  int dosin = flags&S12_FLAG_DOSIN;


  // Accumulators
  double A = 0;
  double B = 0;
  double C = 0;
  for (n_part=0 ; n_part < SD1->N ; n_part++) {
    if ((type != S12_TYPE_ANY) && (SD1->atomtype[n_part] != type))
      continue;
    pos1 = SD1->atompos[n_part];
    pos2 = SD2->atompos[n_part];

    // set dr
    for (d=0 ; d<nDim ; d++) {
      dr[d] =  cords[nDim*pos1 + d] - cords2[nDim*pos2 + d];
      dr[d] -= (floor(dr[d]/lattShape[d] + .5)) *lattShape[d];
    }
    if(print) printf("  dr: %f %f %f (%d)\n", dr[0], dr[1], dr[2], n_part);

    // Loop through k-vectors
    for(n_kvec=0 ; n_kvec < Nk ; n_kvec++) {
      // The regular Fs part.
      double dot=0;
      for (d=0 ; d<nDim ; d++)
	dot += dr[d]*kvecs[nDim*n_kvec+d];
      double x = cos(dot);
/*       double x; */
/*       if (!dosin) */
/* 	x = cos(dot); */
/*       else */
/* 	x = sin(dot); */
      C += x;

      // Now for the prefix part:
      dot = 0;
      for (d=0 ; d<nDim ; d++)
	dot += (cords[nDim*pos1 + d])*q[d];
/*       double y = cos(dot); */
      double y;
      if (!dosin)
	y = cos(dot);
      else
	y = sin(dot);
      A += y*x;
      B += y;

      //SkArrayByAtom[n_part] += x;
      //Skresult[n_kvec] += x;
    }
  }
  *A_ += (double)A/Nk;
  *B_ += (double)B/Nk;
  *C_ += (double)C/Nk;
  return(1);
}



int istructure(struct SimData *SD) {
  // Go through all sites and add particles if we can move.

  int number_of_bonds=0;
  int pos, adjpos;
  int conni;
  for (pos=0 ; pos<SD->lattSize ; pos++) {
    if (SD->lattsite[pos] == S12_EMPTYSITE)
      continue;
    for (conni=0 ; conni<SD->connN[pos] ; conni++) {
      adjpos = SD->conn[pos*SD->connMax + conni];
      if (SD->lattsite[adjpos] != S12_EMPTYSITE)
	number_of_bonds += 1;
    }
  }
  return number_of_bonds;
}


void spinGlass_withOrient(struct SimData *SD0, struct SimData *SD1,
		  int type0, int type1, int *siteCorrelation,int flags);
void spinGlass(struct SimData *SD0, struct SimData *SD1,
		 int type0, int type1,
		 int *siteCorrelation,
		 int flags) {
  int lattSize = SD0->lattSize;
  int pos0, pos1;
  //int hasOrient=0;
  if (SD0->orient != NULL && SD1->orient != NULL) {
    //hasOrient = 1;
    return spinGlass_withOrient(SD0,SD1, type0,type1, siteCorrelation, flags);
  }
  for(pos0=0; pos0<lattSize; pos0++){
    // is the i0 atom correct?
    if (SD0->lattsite[pos0] == S12_EMPTYSITE)
      continue;
    if ((type0!=S12_TYPE_ANY) && (type0!=SD0->atomtype[SD0->lattsite[pos0]]))
      continue;
    for(pos1=0; pos1<lattSize; pos1++) {
/*       if (pos0 != pos1) */
/* 	continue; */
      if (SD1->lattsite[pos1] == S12_EMPTYSITE)
	continue;
      if ((type1!=S12_TYPE_ANY) && (type1!=SD1->atomtype[SD1->lattsite[pos1]]))
	continue;
      //if (hasOrient && (SD0->orient[pos0] != SD1->orient[pos1]))
      //  continue;
      siteCorrelation[pos0*lattSize+pos1] += 1;
    }
  }
}
double spinGlass_sumArray(int *data, int size, int nPoints, double mean,
			  double prefactor) {
  int i;
  double result=0;
  double tmp;
  //int sum=0;
  for (i=0 ; i<size ; i++) {
    //sum += data[i];
    //printf("%d %d %d\n", i, data[i], sum);
    tmp = (((double)data[i])/nPoints) - mean;
    tmp *= tmp;
    result += tmp;
  }
  //printf("sum: %d\n", sum);
  return (prefactor * result);
}
void spinGlass_withOrient(struct SimData *SD0, struct SimData *SD1,
		 int type0, int type1,
		 int *siteCorrelation,
		 int flags) {
  int lattSize = SD0->lattSize;
  int pos0, pos1;
  int nOrient = SD0->connMax;
  for(pos0=0; pos0<lattSize; pos0++){
    // is the i0 atom correct?
    if (SD0->lattsite[pos0] == S12_EMPTYSITE)
      continue;
    if ((type0!=S12_TYPE_ANY) && (type0!=SD0->atomtype[SD0->lattsite[pos0]]))
      continue;
    for(pos1=0; pos1<lattSize; pos1++) {
      if (SD1->lattsite[pos1] == S12_EMPTYSITE)
	continue;
      if ((type1!=S12_TYPE_ANY) && (type1!=SD1->atomtype[SD1->lattsite[pos1]]))
	continue;
      //if (hasOrient && (SD0->orient[pos0] != SD1->orient[pos1]))
      //  continue;
      int index0 = pos0*nOrient+SD0->orient[pos0];
      int index1 = pos1*nOrient+SD1->orient[pos1];
      siteCorrelation[lattSize*nOrient*index0+index1] += 1;
    }
  }
}

void fourpointDensity(struct SimData *SD0, struct SimData *SD1,
		      int type,
		      int *siteCorr4, int *siteCorr2,
		      int flags) {
  int lattSize = SD0->lattSize;
  int pos0, pos1;
  for(pos0=0; pos0<lattSize; pos0++){
    // is the i0 atom correct?
    if (SD0->lattsite[pos0] == S12_EMPTYSITE ||
	SD1->lattsite[pos0] == S12_EMPTYSITE)
      continue;
    if ((type!=S12_TYPE_ANY) && (
				 type!=SD0->atomtype[SD0->lattsite[pos0]] ||
				 type!=SD1->atomtype[SD1->lattsite[pos0]]))
      continue;
    siteCorr2[pos0] += 1;
    for(pos1=0; pos1<lattSize; pos1++) {
      if (SD0->lattsite[pos1] == S12_EMPTYSITE ||
	  SD1->lattsite[pos1] == S12_EMPTYSITE)
	continue;
      if ((type!=S12_TYPE_ANY) && (
				   type!=SD0->atomtype[SD0->lattsite[pos1]] ||
				   type!=SD1->atomtype[SD1->lattsite[pos1]]))
	continue;
      siteCorr4[pos0*lattSize+pos1] += 1;
    }
  }
}

int _Q(struct SimData *SD0, struct SimData *SD1,
       int type, int flags) {
  int lattSize = SD0->lattSize;
  int pos0, pos1;
  int Q=0;
  for(pos0=0; pos0<lattSize; pos0++){
    // is the i0 atom correct?
    if (SD0->lattsite[pos0] == S12_EMPTYSITE)
      continue;
    if ((type!=S12_TYPE_ANY) && type!=SD0->atomtype[SD0->lattsite[pos0]])
      continue;
    for(pos1=0; pos1<lattSize; pos1++) {
      if (SD1->lattsite[pos1] == S12_EMPTYSITE)
	continue;
      if ((type!=S12_TYPE_ANY) && type!=SD1->atomtype[SD1->lattsite[pos1]])
	continue;
      Q += 1;
    }
  }
  return (Q);
}
int Q(struct SimData *SD0, struct SimData *SD1,
       int type, int flags) {
  int lattSize = SD0->lattSize;
  int pos0;
  int Q=0;
  for(pos0=0; pos0<lattSize; pos0++){
    // is the i0 atom correct?
    if (SD0->lattsite[pos0] == S12_EMPTYSITE ||
	SD1->lattsite[pos0] == S12_EMPTYSITE)
      continue;
    if ((type!=S12_TYPE_ANY) && (type!=SD0->atomtype[SD0->lattsite[pos0]] ||
				 type!=SD1->atomtype[SD1->lattsite[pos0]]))
      continue;
    Q += 1;
  }
  return (Q);
}







void ctest(struct SimData *SD) {
  int i;

  printf("N: %d\n", SD->N);
  printf("printing atomtype\n");
  //longX *atomtype = SD->atomtype;
  for (i=0 ; i<SD->lattSize ; i++) {
    //printf("%d ", atomtype[i]);
    printf("%d ", SD->atomtype[i]);
  } printf("\n");

  printf("printing lattsite\n");
  //longX *lattsite = SD->lattsite;
  for (i=0 ; i<SD->lattSize ; i++) {
    //printf("%d ", lattsite[i]);
    printf("%d ", SD->lattsite[i]);

  } printf("\n");
}

int ctest_array(int num, struct SimData **SD_array) {
  int i;
  for (i=0; i<num; i++) {
    struct SimData *SD;
    SD = SD_array[i];
    printf("This is System struct #%d, and has %d atoms\n", i, SD->N);
  }
  return (num);
}
