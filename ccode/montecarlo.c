/* Richard Darst, March 2010 */


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
  // Skip if site is frozen
  if ((SD->flags & S12_FLAG_FROZEN) && SD->frozen[pos])
    return(0);
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
  // Skip if site is frozen
  if ((SD->flags & S12_FLAG_FROZEN) && SD->frozen[pos])
    return(0);
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
    int frozenEnabled = SD->flags & S12_FLAG_FROZEN;
    // Skip if site is frozen
    if (frozenEnabled && SD->frozen[pos])
      return(0);

    // Find an arbitrary connection:
    int i_conn = SD->connN[pos] * genrand_real2();
    //int *connLocal = SD->conn + pos*SD->connMax;
    int newPos = SD->conn[ pos*SD->connMax + i_conn];
    // Skip if adjecent site is frozen
    if (frozenEnabled && SD->frozen[newPos])
      return(0);

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
  // Are there features enabled which if we don't know about them,
  // should cause an error?
  if (SD->flags & S12_FLAG_INCOMPAT & ~S12_FLAG_FROZEN ) {
    printf("Incompatible features seen: %d\n", SD->flags);
    exit(216);
  }

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

