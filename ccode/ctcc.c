/* Richard Darst, May 2009 */

inline int cycleCTCC_translate(struct SimData *SD) {
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

    int oldOrient = -1;
    double Eold, Enew;
    int newPos = -1;
    if (i_conn == SD->orient[pos]) {
      // translation.
      newPos = SD->conn[ pos*SD->connMax + i_conn];
      if (debug) printf("%d %d %d\n", pos, i_conn, newPos);
      if (SD->lattsite[newPos] != S12_EMPTYSITE) {
	// can't move to an adjecent occupied site, reject move.
	if(debug) printf("can't move to adjecent occpuied site\n");
	return(0);
      }
      Eold = energy_pos(SD, pos) + energy_pos(SD, newPos);
      moveParticle(SD, pos, newPos);
      Enew = energy_pos(SD, newPos) + energy_pos(SD, pos);
    }
    else {
      // Internal vibration.
      oldOrient = SD->orient[pos];
      Eold = energy_pos(SD, pos);
      SD->orient[pos] = i_conn;
      Enew = energy_pos(SD, pos);
    }

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
      if (oldOrient == -1)
	moveParticle(SD, newPos, pos);
      else
	SD->orient[pos] = oldOrient;
    }
    else {
      if (debug) printf("accepting move\n");
      if (SD->persist != NULL) { // Update persistence function array if there
        if (oldOrient == -1) {
	  SD->persist[pos] = 1;
	  SD->persist[newPos] = 1;
	}
      }
      return(1);  // Return 1, since we accepted one move
    }
    if(debug) printf("Eold: %.3f Enew: %.3f\n", Eold, Enew);
    return(0);    // Return zero, representing accepting no move
}



int cycleCTCC(struct SimData *SD, double n) {

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
      naccept += cycleCTCC_translate(SD);
    }
  }
  return(naccept);
}







void EddCTCC_init(struct SimData *SD) {
  int pos;
  for (pos=0 ; pos<SD->lattSize ; pos++) {
    if (debugedd)
      printf("pos: %d\n", pos);
    if (SD->lattsite[pos] != S12_EMPTYSITE)
      EddCTCC_updateLatPos(SD, pos);
  }
}

inline void EddCTCC_updateLatPos(struct SimData *SD, int pos) {
  // iterate through all connections, see if each is correctly satisfied...
  //
  int conni;
  if (debugedd)
    printf("EddCTCC_updateLatPos: pos:%d\n", pos);
  int orient = SD->orient[pos];

  for(conni=0 ; conni < SD->connN[pos] ; conni++) {
    // Is this the move where we shift positions?
    int adjpos = SD->conn[SD->connMax*pos + conni];
    if (conni == orient) {
      if (debugedd)
	printf("  adjpos: %d\n", adjpos);
      if (SD->lattsite[adjpos] != S12_EMPTYSITE) {
	ensureNotInMLL(SD, pos, conni);
	continue;
      }
      int allowed_move_to_adjecent = 1;
      int connii;
      for(connii=0 ; connii < SD->connN[adjpos] ; connii++) {
	int adj2pos = SD->conn[SD->connMax*adjpos + connii];
	if (adj2pos == pos)
	  continue;
	if (SD->orient[adj2pos] == reverseConnection(SD, connii)) {
	  allowed_move_to_adjecent = 0;
	  break;
	}
      }
      ensureInMLLIf(allowed_move_to_adjecent, SD, pos, conni);
      continue;
    }
    // So this is a orientation shift only.
    if (SD->lattsite[adjpos] == S12_EMPTYSITE) {
      ensureInMLL(SD, pos, conni);
    } else {
      ensureNotInMLL(SD, pos, conni);
    }
  }
}
inline void EddCTCC_updateLatPos2(struct SimData *SD, int pos,
				struct LList *llist) {
  if (LlistLookup(llist, pos))
    // We are already in the list, so we were already processed--
    // don't process it again.
    return;
  LlistAdd(llist, pos);
  EddCTCC_updateLatPos(SD, pos);
}


int EddCTCC_consistencyCheck(struct SimData *SD) {
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
	// there is a particle here.
	moveIndex = pos*connMax + conni;
	if (conni == SD->orient[pos]) {
	  // We are oriented this way.  Can we move?
	  int adjpos = SD->conn[moveIndex];
	  if (SD->lattsite[adjpos] != S12_EMPTYSITE){
	    // not empty adjecent position, impossible to move there.
	    if (SD->MLLr[moveIndex] != -1) {
	      retval += 1;
	      printf("error vkgcoovocu\n");
	    }
	  } else {
	    // We can possibly move here.  Try and see.
	    moveParticle(SD, pos, adjpos);
	    if (energy_pos(SD, pos) + energy_pos(SD, adjpos) == 1/0. ) {
	      // can't move
	      if (SD->MLLr[moveIndex] != -1) {
		retval += 1;
		printf("error kvocrkn\n");
	      }
	    } else {
	      // can move.
	      if (SD->MLLr[moveIndex] == -1) {
		retval += 1;
		printf("error gcknocgokr\n");
	      }
	    }
	    moveParticle(SD, adjpos, pos);
	  }
	  continue;
	}
	// not oriented this way.  Can we rotate to this angle?
	int oldOrient = SD->orient[pos];
	SD->orient[pos] = conni;
	if ((energy_pos(SD, pos) == 1/0.) && (SD->MLLr[moveIndex] != -1)) {
	  retval += 1;
	  printf("error jrotaaor, pos:%d conni:%d\n",
		 pos, conni);
	}
	if (energy_pos(SD, pos) != 1/0. && SD->MLLr[moveIndex] == -1) {
	  retval += 1;
	  printf("error rkotrouwk, pos:%d conni:%d\n",
		 pos, conni);
	}
	SD->orient[pos] = oldOrient;
      }
    }
  }
  return(retval);
}
int EddCTCC_cycle(struct SimData *SD, double n) {
  if (SD->MLLlen == 0) {
    printf("EddCTCC_cycle: error, move list length is zero\n");
    exit(125);
  }
  int connMax = SD->connMax;
  int naccept = 0;
  struct LList llist; llist.n = 0;

  double maxTime = (double) n;
  double time = SD->MLLextraTime;
  if (time == -1.) {
    // pre-move, advance time until the first event.  Otherwise we
    // always end up moving right at time == 0
    double timestep = (SD->N * SD->connMax) / ((double)SD->MLLlen);
    timestep *= -log(genrand_real3());
    time = timestep;
  }

  //printf("eddCTCC_cycle: time:%f maxTime%f\n", time, maxTime);
  while (time < maxTime) {
    //printf("time:%f\n", time);

    int movei = SD->MLL[(int)(SD->MLLlen*genrand_real2())];
    int oldpos = movei / connMax;
    int conni  = movei % connMax;
    int newpos = SD->conn[movei]; // movei = connMax*oldpos + moveConni
    int oldOrient = -1;
    if (debugedd)
      printf("edd move: moving from oldpos:%d, conni:%d\n", oldpos, conni);

    if (conni == SD->orient[oldpos]) {
      // translation.
      if (debugedd)
	if (SD->lattsite[newpos] != S12_EMPTYSITE) {
	  // can't move to an adjecent occupied site, reject move.
	  printf("error vcoknorourao\n");
      }
      moveParticle(SD, oldpos, newpos);
      if (SD->persist != NULL) {  // Update persistence function array if there
	SD->persist[oldpos] = 1;
	SD->persist[newpos] = 1;
      }
    }
    else {
      // Internal vibration.
      oldOrient = SD->orient[oldpos];
      SD->orient[oldpos] = conni;
    }

    //llist.n = 0;
    //LlistAdd(&llist, oldpos); //oldpos moves are rmvd below, in the first loop.
    int pos_conni;
    if ( oldOrient == -1 ) {
      // a translation move.
      // For all connections of the old site.
      for(pos_conni=0 ; pos_conni < SD->connN[oldpos] ; pos_conni++) {
	// Remove option of oldpos rotating to this connection.
	ensureNotInMLL(SD, oldpos, pos_conni);
	int adjpos = SD->conn[connMax*oldpos+pos_conni];
	// The neighbors of the old site can possibly rotate to face oldpos now
	if (SD->lattsite[adjpos] != S12_EMPTYSITE) {
	  EddCTCC_updateLatPos(SD, adjpos);
	}
      }
      // Now update the new site
      EddCTCC_updateLatPos(SD, newpos);
      // Update all neighbors of the new position
      for(pos_conni=0 ; pos_conni < SD->connN[newpos] ; pos_conni++) {
	int adjpos = SD->conn[connMax*newpos+pos_conni];
	if (SD->lattsite[adjpos] != S12_EMPTYSITE) {
	  EddCTCC_updateLatPos(SD, adjpos);
	}
      }
    }
    else {
      // a vibration move.
      //removeFromMLL(SD, oldpos, conni);
      //addToMLL(SD, oldpos, oldOrient);
      EddCTCC_updateLatPos(SD, oldpos);
      // Update second neighbors of adjpos from oldorient
      int adjpos = SD->conn[connMax*oldpos+conni];
      for(pos_conni=0 ; pos_conni < SD->connN[adjpos] ; pos_conni++) {
	int adj2pos = SD->conn[connMax*adjpos+pos_conni];
	if (SD->lattsite[adj2pos] != S12_EMPTYSITE) {
	  EddCTCC_updateLatPos(SD, adj2pos);
	}
      }
      // Update second neighbors of adjpos from neworient
      adjpos = SD->conn[connMax*oldpos+oldOrient];
      for(pos_conni=0 ; pos_conni < SD->connN[adjpos] ; pos_conni++) {
	int adj2pos = SD->conn[connMax*adjpos+pos_conni];
	if (SD->lattsite[adj2pos] != S12_EMPTYSITE) {
	  EddCTCC_updateLatPos(SD, adj2pos);
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
