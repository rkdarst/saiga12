/* Richard Darst, March 2010 */


inline int cycleSpinMC_flip(struct SimData *SD) {
    // otherwise, do a regular move (this should be the most common
    // case and thus inlined)

    int pos;
    pos = (int)(SD->lattSize * genrand_real2());
    int frozenEnabled = SD->flags & S12_FLAG_FROZEN;
    // Skip if site is frozen
    if (frozenEnabled && SD->frozen[pos])
      return(0);

    int state = SD->lattsite[pos] != S12_EMPTYSITE;
    int newstate = ! state;
    if (debug) printf("%d %d %d\n", pos, state, newstate);

    int accept;
    double Eold = energy_pos(SD, pos);
    if (state) {
      delParticle(SD, pos);
    } else {
      addParticle(SD, pos, SD->inserttype);
    }
    double Enew = energy_pos(SD, pos);
    //printf("%f %f %f\n", Eold, Enew, energy(SD));

    if (Enew == 1./0.) {
      // always reject moves producing infinite energy
      if (debug) printf("illegal move, energy becomes infinite\n");
      accept = 0;
    }
    else if (Enew <= Eold) {
      //printf("Accepting, energy decreasing\n");
      // always accept energy decreasing moves
      accept = 1; }
    else {
      // accept increasing energy moves with metropolis criteria
      double x;
      x = exp(SD->beta*(Eold-Enew));
      double ran = genrand_real2();
      //printf("Testing: rand:%f x:%f arg:%f", ran, x, SD->beta*(Eold-Enew));
      if (ran < x) {
	//printf(" ... accepted\n");
        accept = 1; }
      else {
	//printf(" ... rejected\n");
	accept = 0; }
    }

    if (accept == 0) {
      // Restore the particle location if it wasn't accepted.
      if (state) {
	addParticle(SD, pos, SD->inserttype);
      } else {
	delParticle(SD, pos);
      }
    }
    else {
      if (debug) printf("accepting move\n");
      if (SD->persist != NULL) { // Update persistence function array if there
        SD->persist[pos] = 1;
      }
      return(1);  // Return 1, since we accepted one move
    }
    if(debug) printf("Eold: %.3f Enew: %.3f\n", Eold, Enew);
    return(0);    // Return zero, representing accepting no move
}

int cycleSpinMC(struct SimData *SD, double n) {

  int i_trial;
  int naccept = 0;
  // Are there features enabled which if we don't know about them,
  // should cause an error?
  if (SD->flags & S12_FLAG_INCOMPAT & ~S12_FLAG_FROZEN ) {
    printf("Incompatible features seen: %d\n", SD->flags);
    exit(216);
  }

  for (i_trial=0 ; i_trial<n ; i_trial++) {
    naccept += cycleSpinMC_flip(SD);
  }
  return(naccept);
}

