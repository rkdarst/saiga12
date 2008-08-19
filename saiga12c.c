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
#define S12_TYPE_ANY (-1)

int debug = 0;
int errorcheck = 0;  // print errors if inconsistent
int debugedd = 0;    // print out info always

struct SimData {
  double beta;
  int N;
  int *ntype;
  //int NMax;  // lattSize is NMax
  double hardness;

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
  double MLLextraTime;

  int lattSize;
  int *lattsite;
  int *conn;
  int *connN;
  int connMax;
  int *nneighbors;
  int *atomtype;
  int *atompos;

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


#define neighlist 
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


#ifdef neighlist
inline void addParticle(struct SimData *SD, int pos, int type) {
  if (errorcheck) if (SD->lattsite[pos] != S12_EMPTYSITE) exit(61);
  
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
  if (errorcheck) if (SD->lattsite[pos] == S12_EMPTYSITE) exit(62);
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
inline void moveParticle(struct SimData *SD, int oldpos, int newpos) {
  if (errorcheck) if (SD->lattsite[newpos] != S12_EMPTYSITE) exit(63);
  if (errorcheck) if (SD->lattsite[oldpos] == S12_EMPTYSITE) exit(64);
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
#else

        inline void addParticle(struct SimData *SD, int pos, int type) {
          SD->lattsite[pos] = type;
        }
        inline void delParticle(struct SimData *SD, int pos) {
          SD->lattsite[pos] = S12_EMPTYSITE;
        }

        //#define addParticle(x, y, z) (x)->lattsite[y] = (z);
        // #define delParticle(x, y) (x)->lattsite[y] = S12_EMPTYSITE;

#endif



void loadStateFromSave(struct SimData *SD) {
  int n, pos;
  for (n=0 ; n<SD->N ; n++) {
    SD->lattsite[SD->atompos[n]] = n;
    SD->ntype[SD->atomtype[n]] += 1;
  }
  for (pos=0 ; pos<SD->lattSize ; pos++) {
    SD->nneighbors[pos] = neighbors_pos(SD, pos) ;
  }
}



inline double energy_pos(struct SimData *SD, int pos) {
  /*  Energy of a particle at one particular position.
   */
  if (SD->lattsite[pos] == S12_EMPTYSITE)
     return 0;
  int type = atomType(SD, pos);
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
   * that position, and all neighboring positions.  This takes into
   * account "induced" interactions, where one particle doesn't
   * violate its own density constraint, but violates the density
   * constraint one of its neighbors.
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
  /* Total energy of the system.
   */
  int pos;
  double energy=0;
  for (pos=0 ; pos<SD->lattSize ; pos++) {
    energy += energy_pos(SD, pos);
  }
  return energy;
}





int grandcanonical_add(struct SimData *SD, int pos) {
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
  if (errorcheck) if (SD->lattsite[pos] != S12_EMPTYSITE) exit(57);
  
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
  /* Grand-canonical move where we try to remove a particle.  Called
   * from cycle().
   */
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
  if (errorcheck) if (SD->lattsite[pos] == S12_EMPTYSITE) exit(56);

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
  delParticle(SD, pos);  // this WILL result in particle nums being shifted.
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
  int naccept = 0;

  for (i_trial=0 ; i_trial<n ; i_trial++) {

    //int N = SD->N;

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



    // otherwise, do a regular move (this should be the most common
    // case and thus inlined)

    // Find a lattice site with a particle:
    if (SD->N == 0) {
      printf("we are out of particles!\n");
      exit(165);
    }
    int pos;
    do {
      pos = SD->lattSize * genrand_real2();
    } while (SD->lattsite[pos] == S12_EMPTYSITE);
    //pos = SD->atompos[ (int)(SD->N * genrand_real2()) ];


    // Find an arbitrary connection:
    int i_conn = SD->connN[pos] * genrand_real2();
    //int *connLocal = SD->conn + pos*SD->connMax;
    int newPos = SD->conn[ pos*SD->connMax + i_conn];
    
    if (debug) printf("%d %d %d\n", pos, i_conn, newPos);
    
    if (SD->lattsite[newPos] != S12_EMPTYSITE) {
      // can't move to an adjecent occupied site, reject move.
      if(debug) printf("can't move to adjecent occpuied site\n");
      continue;
    }
    
    int accept;
    double Eold = energy_posNeighborhood(SD, pos) +
                  energy_posNeighborhood(SD, newPos);
    //int atomtype = atomType(SD, pos);
    //delParticle(SD, pos);
    //addParticle(SD, newPos, atomtype);
    moveParticle(SD, pos, newPos);
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
      // Restore the particle location if it wasn't accepted.
      //atomtype = atomType(SD, newPos);
      //delParticle(SD, newPos);
      //addParticle(SD, pos, atomtype);
      moveParticle(SD, newPos, pos);
    }
    else {
      if (debug) printf("accepting move\n");
      naccept += 1;
    }

    if(debug) printf("Eold: %.3f Enew: %.3f\n", Eold, Enew);
  }
  return(naccept);
}




double calc_structfact(struct SimData *SD1, struct SimData *SD2,
		       double *kvecs, int Nk, int type,
		       int *cords, 
		       double *lattShape, int nDim,
		       double *Skresult) {
  int n1, n2, nk;
  double totalsum = 0;
  int pos1, pos2;
  //double r1[3];
  double dr[5];  // max number of dims
  int print=0;

  for (n1=0 ; n1 < SD1->N ; n1++) {
    //printf("%X %X\n", SD1, SD2);
    if ((type != S12_TYPE_ANY) && (SD1->atomtype[n1] != type))
      continue;
    pos1 = SD1->atompos[n1];
    //for (n2=n1+1 ; n2 < SD->N ; n2++) {
    n2 = n1;
    //if ((type != S12_TYPE_ANY) && (SD->atomtype[n2] != type))
    //  continue;
    pos2 = SD2->atompos[n2];
    if (pos1 == pos2) {
      totalsum += Nk;
      for(nk=0 ; nk < Nk ; nk++)
	Skresult[nk] += 1.;
      continue;
    }
    //if (pos1 == pos2)
    //  print = 0;    else print=1;
    if(print) printf("  n1:%d  n2:%d  pos1:%d  pos2:%d\n", n1, n2, pos1, pos2);
    if(print) printf("  cords:  %d %d %d   %d %d %d\n", 
		     cords[3*pos1+0], cords[3*pos1+1], cords[3*pos1+2],
		     cords[3*pos2+0], cords[3*pos2+1], cords[3*pos2+2]);
/*     dr[0] = cords[3*pos1 + 0] - cords[3*pos2 + 0]; */
/*     dr[1] = cords[3*pos1 + 1] - cords[3*pos2 + 1]; */
/*     dr[2] = cords[3*pos1 + 2] - cords[3*pos2 + 2]; */
/*     if(print) printf("  dr: %f %f %f\n", dr[0], dr[1], dr[2]); */
/*     dr[0] -= (floor(dr[0]/lattShape[0] + .5)) *lattShape[0]; */
/*     dr[1] -= (floor(dr[1]/lattShape[1] + .5)) *lattShape[1]; */
/*     dr[2] -= (floor(dr[2]/lattShape[2] + .5)) *lattShape[2]; */

    if(print) printf("  dr: %f %f %f\n", dr[0], dr[1], dr[2]);
    int d;
    for (d=0 ; d<nDim ; d++) {
      //dr[d] = cords[3*pos1 + d] - cords[3*pos2 + d];
      //dr[d] -= (floor(dr[d]/lattShape[d] + .5)) *lattShape[d];
      dr[d] =  cords[nDim*pos1 + d] - cords[nDim*pos2 + d];
      dr[d] -= (floor(dr[d]/lattShape[d] + .5)) *lattShape[d];
    }
    
	
    for(nk=0 ; nk < Nk ; nk++) {
/*       double dot = dr[0]*kvecs[3*nk + 0] + */
/* 	           dr[1]*kvecs[3*nk + 1] + */
/* 	           dr[2]*kvecs[3*nk + 2] ; */
      double dot=0;
      for (d=0 ; d<nDim ; d++) 
	//dot += dr[d]*kvecs[3*nk+d];
	dot += dr[d]*kvecs[nDim*nk+d];
      
      double x = cos(dot);
      if(print) printf("  dot:%f  x:%f   [%f %f %f]\n", dot, x,
		       kvecs[3*nk+0], kvecs[3*nk+1], kvecs[3*nk+2]);
      totalsum += x;
      Skresult[nk] += x;
    }
    //}
  }
  //printf("totalsum: %f ", totalsum);
  return(totalsum);
}





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
inline void updateMLLatPos(struct SimData *SD, int pos) {
  // iterate through all connections, see if each is correctly satisfied...
  // 
  int conni;
  if (debugedd)
    printf("update: pos:%d\n", pos);

  // added in debugging
  //if (SD->lattsite[pos] == S12_EMPTYSITE) {
  //  for(conni=0 ; conni < SD->connN[pos] ; conni++) {
  //    if (SD->MLLr[pos*SD->connMax + conni] != -1) {
  //	removeFromMLL(SD, pos, conni);
  //    }
  //  }
  //  return;
  //}

  for(conni=0 ; conni < SD->connN[pos] ; conni++) {
    int adjpos = SD->conn[SD->connMax*pos + conni];
    if (debugedd)
      printf("adjpos: %d\n", adjpos);
    int isAllowedMove = 1;
    if (SD->lattsite[adjpos] != S12_EMPTYSITE)
      isAllowedMove = 0;
    else {
      // we know it is empty, try moving and see if energy becomes inf.
      moveParticle(SD, pos, adjpos);
      if (energy_posNeighborhood(SD, adjpos) == 1/0.)
	isAllowedMove = 0;
      moveParticle(SD, adjpos, pos);
    }
    if (isAllowedMove) {
      // move is allowed, be sure it is in the lists.
      //printf("x: %d\n", SD->MLLr[pos*SD->connMax + conni]);
      if (SD->MLLr[pos*SD->connMax + conni] == -1)
	addToMLL(SD, pos, conni);
    }
    else {
      // move isn't allowed, remove it if it is in there
      if (SD->MLLr[pos*SD->connMax + conni] != -1)
	removeFromMLL(SD, pos, conni);
    }
  }
}
inline void updateMLLatPos2(struct SimData *SD, int pos, struct LList *llist) {
  if (LlistLookup(llist, pos))
    // We are already in the list, so we were already processed--
    // don't process it again.
    return;
  LlistAdd(llist, pos);
  updateMLLatPos(SD, pos);
}



void initMLL(struct SimData *SD) {
  int pos;
  for (pos=0 ; pos<SD->lattSize ; pos++) {
    if (debugedd)
      printf("pos: %d\n", pos);
    if (SD->lattsite[pos] != S12_EMPTYSITE)
      updateMLLatPos(SD, pos);
  }
}
int MLLConsistencyCheck(struct SimData *SD) {
  // first check that all things in MLLr map to the right thing in MLL
  int MLLlocation;
  int moveIndex;
  int connMax = SD->connMax;
  int retval=0;
  for (moveIndex=0 ; moveIndex<SD->lattSize * connMax ; moveIndex++) {
    if (SD->MLLr[moveIndex] != -1) {
      // if it exists in MLLr, it should be at that point in MLL
      if (moveIndex != SD->MLL[SD->MLLr[moveIndex]] ) {
	retval += 1;
	printf("error orjlhc\n");
      }
    }
  }
  // Now check that all things in the MLL map to the right thing in
  // the MLLr
  for (MLLlocation=0 ; MLLlocation<(SD->lattSize*connMax) ; MLLlocation++) {
    if (MLLlocation < SD->MLLlen) {
      // if it's less than the list length, then it should be look-up able.
      if (MLLlocation != SD->MLLr[SD->MLL[MLLlocation]]) {
	retval += 1;
	printf("error mcaockr\n");
      }
    }
    else {
      // all these greater ones should be blank
      if (SD->MLL[MLLlocation] != -1) {
	retval += 1;
	printf("error rcaohantohk\n");
      }
    }
    
  }
  // Now look and see if everything in MLLr is there if it needs to
  // be... 
  int pos, conni;
  for (pos=0 ; pos<SD->lattSize ; pos++) {
    if (SD->lattsite[pos] == S12_EMPTYSITE) {
      // be sure that it is not in any of the lookups.
      for (conni=0 ; conni<SD->connMax ; conni++) {
	if (SD->MLLr[pos*connMax + conni] != -1) {
	  retval += 1;
	  printf("error aroork\n");
	}
      }
    } else {
      // So we do have a particle here.  Be sure that it is correct...
      // basically reproduce the logic of the update function.
      for (conni=0 ; conni<SD->connMax ; conni++) {
	if (conni >= SD->connN[pos]) {
	  // if this is above our number of connections, it must be empty
	  if (SD->MLLr[pos*connMax + conni] != -1) {
	    retval += 1;
	    printf("error pvwho\n");
	  }
	  continue;
	}
	// this is a real connection.
	moveIndex = pos*connMax + conni;
	int adjpos = SD->conn[moveIndex];
	if (SD->lattsite[adjpos] != S12_EMPTYSITE) {
	  // there is a neighboring particle, this move is not allowed.
	  if (SD->MLLr[pos*connMax + conni] != -1) {
	    retval += 1;
	    printf("error jrotaaor, pos:%d conni:%d adjpos:%d\n",
		   pos, conni, adjpos);
	  }
	} else {
	  // there is not a neighboring particle, so we have to do a
	  // move test.
	  moveParticle(SD, pos, adjpos);
	  if (energy_posNeighborhood(SD, adjpos) == 1/0.) {
	    // not allowed to move here
	    if (SD->MLLr[moveIndex] != -1) {
	      retval += 1;
	      printf("error rcxaorhcu\n");
	    }
	  }
	  else {
	    // we are allowed to move there...
	    if (SD->MLLr[moveIndex] == -1) {
	      retval += 1;
	      printf("error ycorkon\n");
	      printf("... pos:%d adjpos:%d, conni:%d\n", pos, adjpos, conni);
	    }
	  }
	  moveParticle(SD, adjpos, pos);
	}
      }
    }
  }
  return(retval);
}


int eddCycle(struct SimData *SD, int n) {
  
  int connMax = SD->connMax;
  int naccept = 0;
  struct LList llist; llist.n = 0;

  double maxTime = (double) n;
  double time = SD->MLLextraTime;
  //printf("eddCycle: time:%f maxTime%f\n", time, maxTime);
  while (time < maxTime) {

    int movei = SD->MLL[(int)(SD->MLLlen*genrand_real2())];
    int oldpos = movei / connMax;
    int newpos = SD->conn[movei]; // movei = connMax*oldpos + moveConni
    if (debugedd)
      printf("move: moving from oldpos:%d to newpos:%d\n", oldpos, newpos);
    moveParticle(SD, oldpos, newpos);  // should always be valid, else
				       // prev prob.
    llist.n = 0;
    LlistAdd(&llist, oldpos);  // these are removed below, in the first loop.
    //updateMLLatPos(SD, newpos);
    updateMLLatPos2(SD, newpos, &llist);


    // Now, to update all data structures.
    int conni;
    // OLDpos neighbors
    for(conni=0 ; conni < SD->connN[oldpos] ; conni++) {
      // this test sees if we used to be able to move somewhere.  If we
      // used to be able to move, we have to remove it now.
      if (SD->MLLr[oldpos*connMax + conni] != -1) {
        removeFromMLL(SD, oldpos, conni);
      }
      // Now handle adjecent particles, which can now move to the center
      // spot.
      int adjpos = SD->conn[oldpos*connMax + conni];
      if (SD->lattsite[adjpos] != S12_EMPTYSITE) {
        // our adjecent position is occupied, so it could move to oldpos
        // now (maybe).  Update it.
        //updateMLLatPos(SD, adjpos);
        updateMLLatPos2(SD, adjpos, &llist);
	//printf(" . updating neighbor: adjpos:%d\n", adjpos);
	// So neighbor is occpied, so it can be relayed through a
	// empty 2nd neighbor to a 3rd neighbor that has a particle in
	// it.
	int connii;
	for(connii=0 ; connii < SD->connN[adjpos] ; connii++) {
	  int adj2pos = SD->conn[adjpos*connMax + connii];
	  if (SD->lattsite[adj2pos] == S12_EMPTYSITE) {
	    //printf(" . examining empty 2nd neighbor: adj2pos:%d\n", adj2pos);
	    int conniii;
	    for(conniii=0 ; conniii < SD->connN[adj2pos] ; conniii++) {
	      int adj3pos = SD->conn[adj2pos*connMax + conniii];
	      if (SD->lattsite[adj3pos] != S12_EMPTYSITE) {
		//printf(" . updating 3rd neighbor: adj3pos:%d\n", adj3pos);
		//updateMLLatPos(SD, adj3pos);
		updateMLLatPos2(SD, adj3pos, &llist);
	      }
	    }
	  }
	}
      }
      //if (0) 1;
      else {
        // The adjecent position was not occupied, so there is no
        // particle there to move.  But the second-neighbors from here
        // *could* move to the adjpos...
        int connii;
        for(connii=0 ; connii < SD->connN[adjpos] ; connii++) {
	  int adj2pos = SD->conn[adjpos*connMax + connii];
	  if (SD->lattsite[adj2pos] != S12_EMPTYSITE) {
	    //updateMLLatPos(SD, adj2pos);
	    updateMLLatPos2(SD, adj2pos, &llist);
	    //printf(" . updating 2nd neighbor: adj2pos:%d\n", adj2pos);
	  }
        }
      }
    }
  
    // NEWpos stuff
    //printf(" .=new neighbors:\n");
    // we already updated the new position of the moved particle above.
    // go over all neighbors...
    for(conni=0 ; conni < SD->connN[newpos] ; conni++) {
      // Handle adjecent particles, which can now move to the center
      // spot.
      int adjpos = SD->conn[newpos*connMax + conni];
      if (SD->lattsite[adjpos] != S12_EMPTYSITE) {
        // our adjecent position is occupied, so our move here might
        // restrict it some... like remove the move to where we are.
        //updateMLLatPos(SD, adjpos);
        updateMLLatPos2(SD, adjpos, &llist);
	//printf(" . updating neighbor: adjpos:%d\n", adjpos);

	// So neighbor is occpied, so it can be relayed through a
	// empty 2nd neighbor to a 3rd neighbor that has a particle in
	// it.
	int connii;
	for(connii=0 ; connii < SD->connN[adjpos] ; connii++) {
	  int adj2pos = SD->conn[adjpos*connMax + connii];
	  if (SD->lattsite[adj2pos] == S12_EMPTYSITE) {
	    //printf(" . examining empty 2nd neighbor: adj2pos:%d\n", adj2pos);
	    int conniii;
	    for(conniii=0 ; conniii < SD->connN[adj2pos] ; conniii++) {
	      int adj3pos = SD->conn[adj2pos*connMax + conniii];
	      if (SD->lattsite[adj3pos] != S12_EMPTYSITE) {
		//printf(" . updating 3rd neighbor: adj3pos:%d\n", adj3pos);
		//updateMLLatPos(SD, adj3pos);
		updateMLLatPos2(SD, adj3pos, &llist);
	      }
	    }
	  }
	}

      }
      //if (0) 1;
      else {
        // The adjecent position was not occupied, so there is no
        // particle there to move.  But the second-neighbors from here
        // *could* move to the adjpos...
        int connii;
        for(connii=0 ; connii < SD->connN[adjpos] ; connii++) {
  	int adj2pos = SD->conn[adjpos*connMax + connii];
  	if (SD->lattsite[adj2pos] != S12_EMPTYSITE) {
	  //printf(" . updating 2nd neighbor: adj2pos:%d\n", adj2pos);
  	  //updateMLLatPos(SD, adj2pos);
  	  updateMLLatPos2(SD, adj2pos, &llist);
	  }
        }
      }
    }
    naccept += 1;
    
    // Advance time
    double timestep = (SD->N * SD->connMax) / ((double)SD->MLLlen);
    timestep *= -log(genrand_real3());  // exponential distribution of times.
                                        // genrand_real3()  -> (0, 1)
    time += timestep;
    //printf("interval: %f\n", (SD->N * SD->connMax) / (double)SD->MLLlen);
  }
  // MLLextraTime is a positive, the amount to increment before our next stop.
  SD->MLLextraTime = time - maxTime;
  //printf("eddCycle: final time:%f maxTime:%f extraTime:%f\n", 
  // time, maxTime, SD->MLLextraTime);

  return(naccept);
}





void ctest(struct SimData *SD) {
  int i;

  printf("printing atomtype\n");
  //longX *atomtype = SD->atomtype;
  for (i=0 ; i<24 ; i++) {
    //printf("%d ", atomtype[i]);
    printf("%d ", SD->atomtype[i]);
  } printf("\n");

  printf("printing lattsite\n");
  //longX *lattsite = SD->lattsite;
  for (i=0 ; i<24 ; i++) {
    //printf("%d ", lattsite[i]);
    printf("%d ", SD->lattsite[i]);

  } printf("\n");
}
