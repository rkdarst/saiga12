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
#define S12_ENERGY_BMnotzero (10)
#define S12_ENERGY_BMimmobile1 (11)
#define S12_ENERGY_BMimmobile1b (12)
#define S12_CYCLE_MC (1)
#define S12_CYCLE_KA (2)
#define S12_CYCLE_FA (3)
#define S12_CYCLE_CTCC (4)
#define S12_FLAG_VIB_ENABLED (1)

int debug = 0;
int errorcheck = 1;  // print errors if inconsistent
int debugedd = 0;    // print out info always

struct SimData {
  double beta;
  int N;
  int *ntype;
  int ntypeMax;
  //int NMax;  // lattSize is NMax
  double hardness;
  int cycleMode;   // see python for definition
  int energyMode;  // see python for definition

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


inline void addParticle(struct SimData *SD, int pos, int type) {
  if (errorcheck) if (SD->lattsite[pos] != S12_EMPTYSITE) { 
      printf("error: inserting atom at not empty site location: %d\n", pos);
      exit(61); }
  if (errorcheck) if (type > SD->ntypeMax) {
    printf("Inserting particle of type greater than ntypeMax\n");
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





int cycleMC_GCadd(struct SimData *SD, int pos) {
  /*  Grand-canonical move where we trial insert a particle.  cycle()
   *  decides when it is time to call this function.
   */
  if (pos == -1) {
    // pick a random position, occupied or not.
    pos = SD->lattSize * genrand_real2();
    if (SD->lattsite[pos] != S12_EMPTYSITE) {
      // already occupied-- new energy is infinite, can't insert.
      return(0);
    }
  }
  if (errorcheck) if (SD->lattsite[pos] != S12_EMPTYSITE) {
      printf("error: GCadd atom at not empty location: %d\n", pos);
      exit(57); }
  
  int inserttype;
  double inserttype_prob, uVTchempotential;
  if (SD->inserttype != S12_EMPTYSITE) { // Single Component
    inserttype = SD->inserttype;
    inserttype_prob = 1.;
    uVTchempotential = SD->uVTchempotential;
  }
  else {  // Multicomponent
    inserttype = getInsertType(SD);
    inserttype_prob = SD->inserttypes_plookup[inserttype];
    uVTchempotential = SD->inserttypes_mulookup[inserttype];
  }

  double Eold = energy_pos(SD, pos);
  addParticle(SD, pos, inserttype);
  double Enew = energy_pos(SD, pos);
  double x = //(SD->lattSize/(double)(SD->N+1)) *   // this for method A only
               exp( SD->beta * (uVTchempotential - Enew + Eold))
               / inserttype_prob;
  int accept;
  if (x >= 1) {
    accept = 1;
  } else {
    double ran = genrand_real2();
    if (ran < x)
      accept = 1;
    else
      accept = 0;
  }
  
  if (accept) {
    //SD->N += 1;
  } else {
    // revert
    delParticle(SD, pos);
  }
  return(0);
}

int cycleMC_GCdel(struct SimData *SD, int pos) {
  /* Grand-canonical move where we try to remove a particle.  Called
   * from cycle().
   */
  if (pos == -1) {
    // pick a random particle.
    if (SD->N == 0) {
      printf("eNo particles remaining (grandcanonical_del), exiting\n");
      exit(45);
    }
/*     do { */
/*       pos = SD->lattSize * genrand_real2(); */
/*     } while (SD->lattsite[pos] == S12_EMPTYSITE); */
    pos = SD->atompos[ (int)(SD->N * genrand_real2()) ];
  }
  if (errorcheck) if (SD->lattsite[pos] == S12_EMPTYSITE) {
      printf("error: GCdel: removing particle from empty lattsite: %d", pos);
      exit(56); }

  int inserttype;
  double inserttype_prob, uVTchempotential;
  if (SD->inserttype != S12_EMPTYSITE) {  //Singlecomponent
    inserttype_prob = 1.;
    uVTchempotential = SD->uVTchempotential;
  }
  else { // Multicomponent
    inserttype = SD->atomtype[SD->lattsite[pos]];
    inserttype_prob = SD->inserttypes_plookup[inserttype];
    uVTchempotential = SD->inserttypes_mulookup[inserttype];
  }

  double Eold = energy_pos(SD, pos);
  int origtype = atomType(SD, pos);
  int oldOrient = -1;
  if (SD->orient != NULL)
    oldOrient = SD->orient[pos];
  delParticle(SD, pos);  // this WILL result in particle nums being shifted.
  double Enew = energy_pos(SD, pos);
  

  double x = //((SD->N+1)/(double)SD->lattSize) *  // this for method A only
               ((inserttype_prob)) *
               exp( - SD->beta * (uVTchempotential - Eold + Enew));
  int accept;
  if (x >= 1) {
    accept = 1;
  } else {
    double ran = genrand_real2();
    if (ran < x)
      accept = 1;
    else
      accept = 0;
  }
  
  if (accept) {
    //SD->N -= 1;
  } else {
    // revert
    addParticle(SD, pos, origtype);
    if (oldOrient != -1)
      SD->orient[pos] = oldOrient;
  }
  return(0);
}





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



int cycleMC(struct SimData *SD, double n);
int cycleKA(struct SimData *SD, double n);
int cycleFA(struct SimData *SD, double n);
int cycleCTCC(struct SimData *SD, double n);
inline int cycleKA_translate(struct SimData *SD);


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
  else {
    printf("Cycle mode not set: %d", SD->cycleMode);
    exit(49);
  }
}


inline void cycleMC_GCmodeAAdd(struct SimData *SD) {
  // try adding a particle
  //printf("grand canonical: trial add\n");
  cycleMC_GCadd(SD, -1);
}
inline void cycleMC_GCmodeADelete(struct SimData *SD) {
  // try removing a particle
  //printf("grand canonical: trial del\n");
  cycleMC_GCdel(SD, -1);
}

inline void cycleMC_GCmodeB(struct SimData *SD) {
  // try adding a particle
  //printf("grand canonical: trial add\n");
  int pos = SD->lattSize * genrand_real2();
  if (SD->lattsite[pos] == S12_EMPTYSITE)
    cycleMC_GCadd(SD, pos);
  else
    cycleMC_GCdel(SD, pos);
}


inline int cycleMC_translate(struct SimData *SD) {
    // otherwise, do a regular move (this should be the most common
    // case and thus inlined)

    // Find a lattice site with a particle:
    if (SD->N == 0) {
      printf("we are out of particles!\n");
      exit(165);
    }
    int pos;
    pos = SD->atompos[ (int)(SD->N * genrand_real2()) ];


    // Find an arbitrary connection:
    int i_conn = SD->connN[pos] * genrand_real2();
    //int *connLocal = SD->conn + pos*SD->connMax;
    int newPos = SD->conn[ pos*SD->connMax + i_conn];
    
    if (debug) printf("%d %d %d\n", pos, i_conn, newPos);
    
    if (SD->lattsite[newPos] != S12_EMPTYSITE) {
      // can't move to an adjecent occupied site, reject move.
      if(debug) printf("can't move to adjecent occpuied site\n");
      return(0);
      //continue;
    }
    
    int accept;
    double Eold = energy_pos(SD, pos) +
                  energy_pos(SD, newPos);
    //int atomtype = atomType(SD, pos);
    //delParticle(SD, pos);
    //addParticle(SD, newPos, atomtype);
    moveParticle(SD, pos, newPos);
    double Enew = energy_pos(SD, newPos) +
                  energy_pos(SD, pos);

    if (Enew == 1./0.) {
      // always reject moves producing infinite energy
      if (debug) printf("illegal move, energy becomes infinite\n");
      accept = 0;
    }
    else if (Enew <= Eold)
      // always accept energy decreasing moves
      accept = 1;
    else {
      // accept increasing energy moves with metropolis criteria
      double x;
      x = exp(SD->beta*(Eold-Enew));
      double ran = genrand_real2();
      if (ran < x)
        accept = 1;
      else
        accept = 0;
    }
  
    if (accept == 0) {
      // Restore the particle location if it wasn't accepted.
      //atomtype = atomType(SD, newPos);
      //delParticle(SD, newPos);
      //addParticle(SD, pos, atomtype);
      moveParticle(SD, newPos, pos);
    }
    else {
      if (debug) printf("accepting move\n");
      if (SD->persist != NULL) { // Update persistence function array if there
        SD->persist[pos] = 1;
        SD->persist[newPos] = 1;
      }
      return(1);  // Return 1, since we accepted one move
    }
    if(debug) printf("Eold: %.3f Enew: %.3f\n", Eold, Enew);
    return(0);    // Return zero, representing accepting no move
}

int cycleMC(struct SimData *SD, double n) {

  int i_trial;
  int naccept = 0;

  for (i_trial=0 ; i_trial<n ; i_trial++) {

    // Decide what kind of trial move we should do;
    double ran = genrand_real2();

    // method A (pick add or del, then pick a spot needed to make that move)
/*     if (ran < 0 /\*SD->cumProbAdd*\/) { */
/*       cycleMC_GCmodeAAdd(SD); */
/*     } */
/*     else if (ran < 0/\*SD->cumProbDel*\/) { */
/*       cycleMC_GCmodeADelete(SD); */
/*     } */
    // end method A

    // method B  (pick a spot, and then decide if you should try add/del)
    if (ran < SD->cumProbDel) {
      cycleMC_GCmodeB(SD);
    }
    // end method B

    else {
      naccept += cycleMC_translate(SD);
    }
  }
  return(naccept);
}




int cycleKA(struct SimData *SD, double n) {

  int i_trial;
  int naccept = 0;

  for (i_trial=0 ; i_trial<n ; i_trial++) {
    naccept += cycleKA_translate(SD);
  }
  return(naccept);
}

inline int cycleKA_translate(struct SimData *SD) {
    // otherwise, do a regular move (this should be the most common
    // case and thus inlined)

    // Find a lattice site with a particle:
    if (SD->N == 0) {
      printf("we are out of particles!\n");
      exit(165);
    }
    int pos;
    pos = SD->atompos[ (int)(SD->N * genrand_real2()) ];
    int atomtype = atomType(SD, pos);
    if (debug) printf("try kob-andersen move: from %d\n", pos);

    // Is our current site too packed to move from?
    int nneighbors = SD->nneighbors[pos];
    if (nneighbors > atomtype) {
      if(debug) printf("    old site has too many neighbors\n");
      return(0);
    }

    // Find an arbitrary connection:
    int i_conn = SD->connN[pos] * genrand_real2();
    int newPos = SD->conn[ pos*SD->connMax + i_conn];

    if (debug) printf("                     ... to %d\n", newPos);
    // can't move to an adjecent occupied site, reject move:
    if (SD->lattsite[newPos] != S12_EMPTYSITE) {
      if(debug) printf("    can't move to adjecent occpuied site\n");
      return(0);
    }

    // Is the new site too packed to move to?
    nneighbors = SD->nneighbors[newPos];
    if (nneighbors > (atomtype+1)) {
      if(debug) printf("    new site has too many neighbors\n");
      return(0);
    }

    moveParticle(SD, pos, newPos);
    if (SD->persist != NULL) { // Update persistence function array if there
      SD->persist[pos] = 1;
      SD->persist[newPos] = 1;
    }

    return(1);  // Return 1, since we accepted one move
}




double calc_structfact(struct SimData *SD1, struct SimData *SD2,
		       double *kvecs, int Nk, int type,
		       double *cords, double *cords2,
		       double *lattShape, int nDim,
		       double *Skresult,
		       double *SkArrayByAtom,
		       int flags) {
  int n1, n2, nk;
  double totalsum = 0;
  int pos1, pos2;
  //double r1[3];
  double dr[5];  // max number of dims
  int print=0;

  for (n1=0 ; n1 < SD1->N ; n1++) {
    //printf("%p %p\n", SD1, SD2);
    if (print) printf("=== n: %d ===\n", n1);
    if ((type != S12_TYPE_ANY) && (SD1->atomtype[n1] != type))
      continue;
    pos1 = SD1->atompos[n1];
    //for (n2=n1+1 ; n2 < SD->N ; n2++) {
    n2 = n1;
    //if ((type != S12_TYPE_ANY) && (SD->atomtype[n2] != type))
    //  continue;
    pos2 = SD2->atompos[n2];
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
  if (errorcheck && SD->MLLr[lookup] != -1) {
    printf("error rhoeacurk %d %d \n", pos, conni);
  }
  SD->MLL[SD->MLLlen] = lookup;
  SD->MLLr[lookup] = SD->MLLlen;
  SD->MLLlen += 1;
}
inline void removeFromMLL(struct SimData *SD, int pos, int conni) {
  int moveIndex = pos * SD->connMax + conni;
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

inline void EddBM_updateLatPos(struct SimData *SD, int pos);
#include "ccode/birolimezard.c"


inline void EddKA_updateLatPos(struct SimData *SD, int pos);
#include "ccode/kobandersen.c"

inline void EddFA_updateLatPos(struct SimData *SD, int pos);
#include "ccode/fredricksonandersen.c"

#include "ccode/ctcc.c"







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
