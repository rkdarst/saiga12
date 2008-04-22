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

struct SimData {
  double beta;
  int N;
  //int NMax;  // use lattSize instead
  double hardness;
  double uVTchempotential;
  int inserttype;  // must be generalized later.
  int widominserttype;  // must be generalized later.
  int inserttypes_n;    // not used
  double *inserttypes_prob; // when multiple insert types, this is cum prob.
  int *inserttypes_type; // and this is corresponding type.
  double *inserttypes_plookup; // lookup from type->prob
  double *inserttypes_mulookup; // lookup from type->mu
  int lattSize;
  //int *partpos;
  int *lattsite;
  int *conn;
  int *connN;
  int connMax;
  int *nneighbors;
  int *atomtype;
  int *atompos;
  int *ntype;

  double cumProbAdd;
  double cumProbDel;
} ;

/* genrand_real1  [0,1]
 * genrand_real2  [0,1)
 * genrand_real3  (0,0)
 */


#define neighlist 


/* Given lattice point latI, how many adjacent atoms does it have?
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


inline int getInsertType_i(struct SimData *SD) {
  double ran = genrand_real2();  /*[0,1)*/
  int i=0;
  while (1) {
    if (SD->inserttypes_prob[i] >= ran)
      return i;
    i++;
  }
}
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

void printLatticeC(struct SimData *SD) {
  int pos;
  for (pos=0 ; pos<SD->lattSize ; pos++)
    printf("%d ", SD->lattsite[pos]);
  printf("\n");
}

inline int atomType(struct SimData *SD, int pos) {
  //return(SD->lattsite[pos]);
  return(SD->atomtype[SD->lattsite[pos]]);
}


#ifdef neighlist
inline void addParticle(struct SimData *SD, int pos, int type) {
  //if (SD->lattsite[pos] != S12_EMPTYSITE) exit(61);
  
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
  SD->N++;
}
inline void delParticle(struct SimData *SD, int pos) {
  //if (SD->lattsite[pos] == S12_EMPTYSITE) exit(62);
  //SD->lattsite[pos] = S12_EMPTYSITE;
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
    SD->lattsite[pos] = S12_EMPTYSITE;
  }
  else {
    // move the N-1 atom to atomnum.
    SD->atomtype[atomnum] = SD->atomtype[SD->N-1];
    SD->atompos[atomnum] = SD->atompos[SD->N-1];
    SD->lattsite[SD->atompos[SD->N-1]] = atomnum;

    SD->atomtype[SD->N-1] = S12_EMPTYSITE;  // reset it
    SD->atompos[SD->N-1] = S12_EMPTYSITE;  // reset it

    SD->lattsite[pos] = S12_EMPTYSITE;
  }
  SD->N--;
}
#else

        inline void addParticle(struct SimData *SD, int pos, int type) {
          SD->lattsite[pos] = type;
        }
        inline void delParticle(struct SimData *SD, int pos) {
          SD->lattsite[pos] = S12_EMPTYSITE;
        }

        #define addParticle(x, y, z) (x)->lattsite[y] = (z);
        #define delParticle(x, y) (x)->lattsite[y] = S12_EMPTYSITE;

#endif




inline double energy_pos(struct SimData *SD, int pos) {
  // lattice site contiains the type of particle, which is max # of
  // neighbors of that particle.
  int type = atomType(SD, pos);
  if (type == S12_EMPTYSITE)
     return 0;
#ifndef neighlist
  int excessneighbors = neighbors_pos(SD, pos) - type;
#else
  int excessneighbors = SD->nneighbors[pos] - type;
#endif
  //printf("excessneighbors: %d\n", excessneighbors);
  if (excessneighbors <= 0) {
    return 0.;
  }
  else {
    return excessneighbors * SD->hardness;
  }
}





inline double energy_posNeighborhood(struct SimData *SD, int pos) {
  /* This is similar to the energy_pos, but will also return energy of
   * that position, and all neighboring positions.
   */
  int i_conn;
  double E = energy_pos(SD, pos);
  for (i_conn=0; i_conn<SD->connN[pos] ; i_conn++) {
    int adjPos = SD->conn[pos*SD->connMax + i_conn];
    E += energy_pos(SD, adjPos);
  }
  return E;
}




double energy(struct SimData *SD) {
  int pos;
  double energy=0;
  for (pos=0 ; pos<SD->lattSize ; pos++) {
    energy += energy_pos(SD, pos);
  }
  return energy;
}





int grandcanonical_add(struct SimData *SD, int pos) {
  if (pos == -1) {
    // pick a random position, occupied or not.
    pos = SD->lattSize * genrand_real2();
    if (SD->lattsite[pos] != S12_EMPTYSITE) {
      // already occupied-- new energy is infinite, can't insert.
      return(0);
    }
  }
  if (SD->lattsite[pos] != S12_EMPTYSITE) exit(57);
  
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

  double Eold = energy_posNeighborhood(SD, pos);
  addParticle(SD, pos, inserttype);
  double Enew = energy_posNeighborhood(SD, pos);
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

int grandcanonical_del(struct SimData *SD, int pos) {
  if (pos == -1) {
    // pick a random particle.
    if (SD->N == 0) {
      printf("eNo particles remaining (grandcanonical_del), exiting\n");
      exit(45);
    }
    do {
      pos = SD->lattSize * genrand_real2();
    } while (SD->lattsite[pos] == S12_EMPTYSITE);
  }
  if (SD->lattsite[pos] == S12_EMPTYSITE) exit(56);

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

  double Eold = energy_posNeighborhood(SD, pos);
  int origtype = atomType(SD, pos);
  delParticle(SD, pos);
  double Enew = energy_posNeighborhood(SD, pos);
  

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
  }
  return(0);
}





double chempotential(struct SimData *SD, int inserttype) {
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
    double Eold = energy_posNeighborhood(SD, pos);
    addParticle(SD, pos, inserttype);
    double Enew = energy_posNeighborhood(SD, pos);
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





int cycle(struct SimData *SD, int n) {

  int i_trial;
  int debug = 0;

  for (i_trial=0 ; i_trial<n ; i_trial++) {

    //int N = SD->N;
    int i_conn;
    int *connLocal;

    // Decide what kind of trial move we should do;
    double ran = genrand_real2();
    // method A (pick add or del, then pick a spot needed to make that move)
/*     if (ran < SD->cumProbAdd) { */
/*       // try adding a particle */
/*       //printf("grand canonical: trial add\n"); */
/*       grandcanonical_add(SD, -1); */
/*       continue; */
/*     } */
/*     else if (ran < SD->cumProbDel) { */
/*       // try removing a particle */
/*       //printf("grand canonical: trial del\n"); */
/*       grandcanonical_del(SD, -1); */
/*       continue; */
/*     } */
    // end method A
    // method B  (pick a spot, and then decide if you should try add/del)
    if (ran < SD->cumProbDel) {
      // try adding a particle
      //printf("grand canonical: trial add\n");
      int pos = SD->lattSize * genrand_real2();
      if (SD->lattsite[pos] == S12_EMPTYSITE)
	grandcanonical_add(SD, pos);
      else
	grandcanonical_del(SD, pos);
      continue;
    }
    // end method B



    // otherwise, do a regular move (this should be the most common case)


    
    //CHECK
    // Find a lattice site with a particle:
    if (SD->N == 0) {
      printf("we are out of particles!\n");
      exit(165);
    }
    int pos;
    do {
      pos = SD->lattSize * genrand_real2();
    } while (SD->lattsite[pos] == S12_EMPTYSITE);
    connLocal = SD->conn + pos*SD->connMax;

    // Find an arbitrary connection:
    i_conn = SD->connN[pos] * genrand_real2();
    int newPos = connLocal[i_conn];
    
    if (debug) printf("%d %d %d\n", pos, i_conn, newPos);
    
    if (SD->lattsite[newPos] != S12_EMPTYSITE) {
      // can't move to an adjecent occupied site.
      if(debug) printf("can't move to adjecent occpuied site\n");
      continue;
    }
    
    int accept;
    double Eold = energy_posNeighborhood(SD, pos) +
                  energy_posNeighborhood(SD, newPos);
/* #ifndef neighlist */
/*     SD->lattsite[newPos] = SD->lattsite[pos]; */
/*     SD->lattsite[pos] = S12_EMPTYSITE; */
/* #else */
    int atomtype = atomType(SD, pos);
    delParticle(SD, pos);
    addParticle(SD, newPos, atomtype);
/* #endif */
    double Enew = energy_posNeighborhood(SD, newPos) +
                  energy_posNeighborhood(SD, pos);

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
      ran = genrand_real2();
      if (ran < x)
        accept = 1;
      else
        accept = 0;
    }
  
    if (accept == 0) {
/* #ifndef neighlist */
/*       SD->lattsite[pos] = SD->lattsite[newPos]; */
/*       SD->lattsite[newPos] = S12_EMPTYSITE; */
/* #else */
      atomtype = atomType(SD, newPos);
      delParticle(SD, newPos);
      addParticle(SD, pos, atomtype);
/* #endif */
    }
    else {
      if (debug) printf("accepting move\n");
    }
    
    
    if(debug) printf("Eold: %.3f Enew: %.3f\n", Eold, Enew);
    
  }
  
  return(0);
}
