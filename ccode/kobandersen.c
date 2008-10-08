void EddKA_init(struct SimData *SD) {
  int pos;
  for (pos=0 ; pos<SD->lattSize ; pos++) {
    if (debugedd)
      printf("pos: %d\n", pos);
    if (SD->lattsite[pos] != S12_EMPTYSITE)
      EddKA_updateLatPos(SD, pos);
  }
}

inline void EddKA_updateLatPos(struct SimData *SD, int pos) {
  // iterate through all connections, see if each is correctly satisfied...
  // 
  int conni;
  if (debugedd)
    printf("EddKA_updateLatPos: pos:%d\n", pos);

  int atomtype = atomType(SD, pos);
  int nneighbors = SD->nneighbors[pos];
  int isAllowedMove_self = 1;
  if (nneighbors > atomtype)
    isAllowedMove_self = 0;


  for(conni=0 ; conni < SD->connN[pos] ; conni++) {
    int adjpos = SD->conn[SD->connMax*pos + conni];
    if (debugedd)
      printf("  adjpos: %d\n", adjpos);
    int isAllowedMove_other = 1;
    if (SD->lattsite[adjpos] != S12_EMPTYSITE)
      isAllowedMove_other = 0;
    else {
      // we know it is empty, is it too packed to move to?
      nneighbors = SD->nneighbors[adjpos]; // nneighbors is changed.
      if (nneighbors > (atomtype+1))
	isAllowedMove_other = 0;
    }
    if (isAllowedMove_self && isAllowedMove_other) {
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
inline void EddKA_updateLatPos2(struct SimData *SD, int pos, 
				struct LList *llist) {
  if (LlistLookup(llist, pos))
    // We are already in the list, so we were already processed--
    // don't process it again.
    return;
  LlistAdd(llist, pos);
  EddKA_updateLatPos(SD, pos);
}


int EddKA_consistencyCheck(struct SimData *SD) {
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

	  // This recalculates nneighbors every time, which is less
          // efficient, but since this is consistencyCheck, it is less
	  // important.
	  int atomtype = atomType(SD, pos);
	  int nneighbors = SD->nneighbors[pos];
	  int isAllowedMove = 1;
	  if (nneighbors > atomtype)
	    isAllowedMove = 0;

	  nneighbors = SD->nneighbors[adjpos];
	  if (nneighbors > (atomtype+1))
	    isAllowedMove = 0;
	  if (isAllowedMove == 0) {
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
	}
      }
    }
  }
  return(retval);
}
int EddKA_cycle(struct SimData *SD, double n) {
  if (SD->MLLlen == 0) {
    printf("EddKA_cycle: error, move list length is zero\n");
    exit(12);
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

  //printf("eddKA_cycle: time:%f maxTime%f\n", time, maxTime);
  while (time < maxTime) {

    int movei = SD->MLL[(int)(SD->MLLlen*genrand_real2())];
    int oldpos = movei / connMax;
    int newpos = SD->conn[movei]; // movei = connMax*oldpos + moveConni
    if (debugedd)
      printf("move: moving from oldpos:%d to newpos:%d\n", oldpos, newpos);
    moveParticle(SD, oldpos, newpos);  // should always be valid, else
				       // prev prob.
    llist.n = 0;
    LlistAdd(&llist, oldpos); //oldpos moves are rmvd below, in the first loop.
    EddKA_updateLatPos2(SD, newpos, &llist);


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
        EddKA_updateLatPos2(SD, adjpos, &llist);
	//printf(" . updating neighbor: adjpos:%d\n", adjpos);
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
	    //EddBM_updateLatPos(SD, adj2pos);
	    EddKA_updateLatPos2(SD, adj2pos, &llist);
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
        //EddBM_updateLatPos(SD, adjpos);
        EddKA_updateLatPos2(SD, adjpos, &llist);
	//printf(" . updating neighbor: adjpos:%d\n", adjpos);
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
  	  //EddBM_updateLatPos(SD, adj2pos);
  	  EddKA_updateLatPos2(SD, adj2pos, &llist);
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
