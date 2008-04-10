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
  int NMax;
  double hardness;
  double chempotentialEx;
  int inserttype;  // must be generalized later.
  int lattSize;
  //int *partpos;
  int *lattsite;
  int *conn;
  int *connN;
  int connMax;

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
    //printf("%d\n", i);
    int neighpos = SD->conn[SD->connMax*pos + i];
    //printf("%d\n", neighpos);
    if (SD->lattsite[neighpos] != S12_EMPTYSITE) {
      neighbors += 1;
    }
  }
  return(neighbors);
}





inline double energy_pos(struct SimData *SD, int pos) {
  // lattice site contiains the type of particle, which is max # of
  // neighbors of that particle.
  if (SD->lattsite[pos] == S12_EMPTYSITE)
     return 0;
  int excessneighbors = neighbors_pos(SD, pos) - SD->lattsite[pos];
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





int grandcanonical_add(struct SimData *SD) {
  /* Insert a particle in a random position
   */
  // pick a random position, occupied or not.
  int pos = SD->lattSize * genrand_real2();
  if (SD->lattsite[pos] != S12_EMPTYSITE) {
    // already occupied-- insert energy is infinite, can't insert.
    return(0);
  }

  double Eold = energy_posNeighborhood(SD, pos);
  SD->lattsite[pos] = SD->inserttype;
  double Enew = energy_posNeighborhood(SD, pos);
  
  int accept;
  double x = SD->beta * //(SD->lattSize/SD->N) *
                        (SD->chempotentialEx - Enew + Eold);
  if (x >= 0) {
    accept = 1;
  } else {
    x = exp(x);
    double ran;
    ran = genrand_real2();
    if (ran < x)
      accept = 1;
    else
      accept = 0;
  }
  
  if (accept) {
    SD->N += 1;
  } else {
    // revert
    SD->lattsite[pos] = S12_EMPTYSITE;
  }
  return(0);
}

int grandcanonical_del(struct SimData *SD) {
  /* Remove a random particle.
   */
  int pos;
  // pick a random particle.
  do {
    pos = SD->lattSize * genrand_real2();
  } while (SD->lattsite[pos] == S12_EMPTYSITE);
  
  double Eold = energy_posNeighborhood(SD, pos);
  int origtype = SD->lattsite[pos];
  SD->lattsite[pos] = S12_EMPTYSITE;
  double Enew = energy_posNeighborhood(SD, pos);
  

  int accept;
  double x = SD->beta * //(SD->N/SD->lattSize) *
                        ((-SD->chempotentialEx) - Enew + Eold);
  if (x >= 0) {
    accept = 1;
  } else {
    x = exp(x);
    double ran;
    ran = genrand_real2();
    if (ran < x)
      accept = 1;
    else
      accept = 0;
  }
  
  if (accept) {
    SD->N -= 1;
  } else {
    // revert
    SD->lattsite[pos] = origtype;
  }
  return(0);
}





double chempotentialEx(struct SimData *SD) {
  // this variable will hold the sum of all exp(-beta deltaU)
  double totalsum=0;
  int inserttype = SD->inserttype;

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
    SD->lattsite[pos] = inserttype;
    double Enew = energy_posNeighborhood(SD, pos);
    SD->lattsite[pos] = S12_EMPTYSITE ;
    //printf("%f %f\n", Eold, Enew);

    totalsum += exp(-SD->beta * (Enew-Eold));
  }
  //printf("%f\n", totalsum);

  double chempotentialEx = - log(totalsum / SD->lattSize) / SD->beta;
  return(chempotentialEx);
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
    if (ran < SD->cumProbAdd) {
      // try adding a particle
      //printf("grand canonical: trial add\n");
      grandcanonical_add(SD);
      continue;
    }
    else if (ran < SD->cumProbDel) {
      // try removing a particle
      //printf("grand canonical: trial del\n");
      grandcanonical_del(SD);
      continue;
    }
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
    double Eold = energy_posNeighborhood(SD, pos);
    SD->lattsite[newPos] = SD->lattsite[pos];
    SD->lattsite[pos] = S12_EMPTYSITE;
    double Enew = energy_posNeighborhood(SD, newPos);

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
      SD->lattsite[pos] = SD->lattsite[newPos];
      SD->lattsite[newPos] = S12_EMPTYSITE;
    }
    else {
      if (debug) printf("accepting move\n");
    }
    
    
    if(debug) printf("Eold: %.3f Enew: %.3f\n", Eold, Enew);
    
  }
  
  return(0);
}
