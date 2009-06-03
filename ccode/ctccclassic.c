/* Richard Darst, May 2009 */

inline int cycleCTCCclassic_translate(struct SimData *SD) {
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
    int nConn = SD->connN[pos];
    int i_conn = (nConn+1) * genrand_real2();
    

    int newPos;
    if (i_conn == nConn) {
      // we are staying in our same location.
      newPos = pos;
    }
    else {
      newPos = SD->conn[ pos*SD->connMax + i_conn];
      if (SD->lattsite[newPos] != S12_EMPTYSITE)
	return(0);
    }
    int oldOrient = SD->orient[pos];
    int newOrient = SD->connN[newPos] * genrand_real2();


    double Eold = energy_pos(SD, pos) +
                  energy_pos(SD, newPos);
    if (pos != newPos)
      moveParticleShort(SD, pos, newPos);
    SD->orient[newPos] = newOrient;
    double Enew = energy_pos(SD, newPos) +
                  energy_pos(SD, pos);

    int accept;
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
      if (pos != newPos)
	moveParticleShort(SD, newPos, pos);
      SD->orient[pos] = oldOrient;
    }
    else {
      if (debug) printf("accepting move\n");
      if (SD->persist != NULL) { // Update persistence function array if there
        if (pos != newPos) {
	  SD->persist[pos] = 1;
	  SD->persist[newPos] = 1;
	}
      }
      return(1);  // Return 1, since we accepted one move
    }
    if(debug) printf("Eold: %.3f Enew: %.3f\n", Eold, Enew);
    return(0);    // Return zero, representing accepting no move
}



int cycleCTCCclassic(struct SimData *SD, double n) {

  int i_trial;
  int naccept = 0;

  for (i_trial=0 ; i_trial<n ; i_trial++) {
    // Decide what kind of trial move we should do;
    double ran = genrand_real2();
    // method B  (pick a spot, and then decide if you should try add/del)
    if (ran < SD->cumProbDel) {
      cycleMC_GCmodeB(SD);
    }
    // end method B

    else {
      naccept += cycleCTCCclassic_translate(SD);
    }
  }
  return(naccept);
}



