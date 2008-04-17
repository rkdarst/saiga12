/* Richard Darst, 2008
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#define CHECK printf(".\n");

#include "SFMT.h"

#define S12_EMPTYSITE 0

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
  int lattSize;
  //int *partpos;
  int *lattsite;
  int *conn;
  int *connN;
  int connMax;
  int *nneighbors;

  double cumProbAdd;
  double cumProbDel;
} ;

/* genrand_real1  [0,1]
 * genrand_real2  [0,1)
 * genrand_real3  (0,0)
 */


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
	return (SD->inserttypes_type[i]);
      i++;
    }
  }
}



#define neighlist 

#ifdef neighlist
inline void addParticle(struct SimData *SD, int pos, int type) {
  //if (SD->lattsite[pos] != S12_EMPTYSITE) exit(61);
  SD->lattsite[pos] = type;
  int i;
  for (i=0 ; i<SD->connN[pos] ; i++) {
    int neighpos = SD->conn[SD->connMax*pos + i];
    SD->nneighbors[neighpos] ++;
  }
}
inline void delParticle(struct SimData *SD, int pos) {
  //if (SD->lattsite[pos] == S12_EMPTYSITE) exit(62);
  SD->lattsite[pos] = S12_EMPTYSITE;
  int i;
  for (i=0 ; i<SD->connN[pos] ; i++) {
    int neighpos = SD->conn[SD->connMax*pos + i];
    SD->nneighbors[neighpos] --;
  }
}
#else

#define addParticle(x, y, z) (x)->lattsite[y] = (z);
#define delParticle(x, y) (x)->lattsite[y] = S12_EMPTYSITE;

#endif




inline double energy_pos(struct SimData *SD, int pos) {
  // lattice site contiains the type of particle, which is max # of
  // neighbors of that particle.
  if (SD->lattsite[pos] == S12_EMPTYSITE)
     return 0;
#ifndef neighlist
  int excessneighbors = neighbors_pos(SD, pos) - SD->lattsite[pos];
#else
  int excessneighbors = SD->nneighbors[pos] - SD->lattsite[pos];
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
  double Eold = energy_posNeighborhood(SD, pos);
  //#ifndef neighlist
  //  SD->lattsite[pos] = SD->inserttype;
  //#else
  addParticle(SD, pos, getInsertType(SD));
  //#endif
  double Enew = energy_posNeighborhood(SD, pos);
  int accept;
  double x = //(SD->lattSize/(double)(SD->N+1)) *   // this for method A only
               exp( SD->beta * ((SD->uVTchempotential) - Enew + Eold));
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
    SD->N += 1;
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
      printf("No particles remaining (grandcanonical_del), exiting\n");
      exit(45);
    }
    do {
      pos = SD->lattSize * genrand_real2();
    } while (SD->lattsite[pos] == S12_EMPTYSITE);
  }
      
  if (SD->lattsite[pos] == S12_EMPTYSITE) exit(56);
  double Eold = energy_posNeighborhood(SD, pos);
  int origtype = SD->lattsite[pos];
  delParticle(SD, pos);
  double Enew = energy_posNeighborhood(SD, pos);
  

  int accept;
  double x = //((SD->N+1)/(double)SD->lattSize) *  // this for method A only
               exp( SD->beta * ((-SD->uVTchempotential) + Eold - Enew));
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
    SD->N -= 1;
  } else {
    // revert
    addParticle(SD, pos, origtype);
  }
  return(0);
}





double chempotential(struct SimData *SD) {
  // this variable will hold the sum of all exp(-beta deltaU)
  double totalsum=0;
  int inserttype = SD->widominserttype;
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
  double chempotential = - log(totalsum / (SD->N+1)) / SD->beta; // B
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
    addParticle(SD, newPos, SD->lattsite[pos]);
    delParticle(SD, pos);
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
      addParticle(SD, pos, SD->lattsite[newPos]);
      delParticle(SD, newPos);
/* #endif */
    }
    else {
      if (debug) printf("accepting move\n");
    }
    
    
    if(debug) printf("Eold: %.3f Enew: %.3f\n", Eold, Enew);
    
  }
  
  return(0);
}
