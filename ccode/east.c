int cycleEast(struct SimData *SD, double n) {

  printf("error: Using East without event driven dynamics does NOT WORK(yet)\n");
  printf("Enable event-driven dynamics (.eddEnable())\n");
  exit(94);

  int i_trial;
  int naccept = 0;

  int N = SD->N;   // initial N value, use SD->N for instantaneous
  double c = 1./(1. + exp(SD->beta));
  //printf("beta:%f c:%f\n", SD->beta, c);
  double wUp = 1-(c);
  double wDown = c;
  //double maxTime = N*n;
  //double time ;

  //for (i_trial=0 ; i_trial<(n*SD->N) ; i_trial++) {
  for (i_trial=0 ; naccept<(N*n) ; i_trial++ ) {

    // otherwise, do a regular move (this should be the most common
    // case and thus inlined)

    // Find a lattice site with a particle:
    if (SD->N == 0) {
      printf("we are out of particles!\n");
      exit(165);
    }
    int pos;
    int atomtype;
    int state;
    while(1) {
      pos = (int)(SD->lattSize * genrand_real2());
      if (SD->lattsite[pos] != S12_EMPTYSITE) {
	// particle present, can we flip it down?
	atomtype = atomType(SD, pos);
	state = 1;
      } else { // trip flipping down to up.
	atomtype = SD->inserttype;
	state = 0;
      }
      if (SD->nneighbors[pos] >= atomtype)
	// we are mobile here-- continue with the move function from this spot.
	break;
    }
    
    double rand = genrand_real2();
    if (state == 1 && rand < wUp) {
      if(debugedd) 
	printf("accepting East move: %d up -> down\n", pos);
      naccept += 1;
      delParticle(SD, pos);
    } else if (state == 0 && rand < wDown) {
      if(debugedd) 
	printf("accepting East move: %d down -> up\n", pos);
      naccept += 1;
      addParticle(SD, pos, SD->inserttype);
    } else {
      if(debugedd) 
	printf("rejecting East move East move: %d (was %d)\n", pos, state);
    }
  }

  return(naccept);  // Return 1, since we accepted one move
}




inline int EastNNeighbors(struct SimData *SD, int pos) {
  //return(SD->nneighbors[pos]);  // Old way
  int neighbors=0;
  int i;
  int nconn_2=SD->connN[pos]/2;
  for (i=0 ; i<nconn_2 ; i++) {
    int neighpos = SD->conn[SD->connMax*pos + i];
    if (SD->lattsite[neighpos] != S12_EMPTYSITE)
      neighbors += 1;
  }
  return (neighbors);
}

void EddEast_init(struct SimData *SD) {
  int pos;
  for (pos=0 ; pos<SD->lattSize ; pos++) {
    if (debugedd)
      printf("EddEast_init: at pos: %d\n", pos);
    EddEast_updateLatPos(SD, pos);
  }
}

inline void EddEast_updateLatPos(struct SimData *SD, int pos) {
  // iterate through all connections, see if each is correctly satisfied...
  // 
  int isAllowedMove=0;
  if (debugedd)
    printf("EddEast_updateLatPos: at update: pos:%d\n", pos);


  int state = SD->lattsite[pos] != S12_EMPTYSITE;
  if(debugedd) 
    printf("state is: %d (should be 0 or 1)\n", state);
  if (state == 1) {
    int atomtype = SD->atomtype[SD->lattsite[pos]];
    if(debugedd) 
      printf("pos:%d nneighbors:%d atomtype:%d\n", pos, 
	     EastNNeighbors(SD, pos), atomtype);
    if (EastNNeighbors(SD, pos) >= atomtype)
      isAllowedMove = 1;
  }
  else if (state == 0  && EastNNeighbors(SD, pos) >= SD->inserttype)
    {
      isAllowedMove = 1;
    }
  else
    isAllowedMove = 0;


  if (isAllowedMove) {
    // move is allowed, be sure it is in the lists.
    if (state == 1) {
      // state is down, add to up lists
      if (SD->MLLr[pos] == -1)
	FAaddToMLL(SD, 'u', pos);
    }
    else {
      // state is down, check the down lists
      if (SD->MLLr[pos] == -1)
	FAaddToMLL(SD, 'd', pos);
    }
  }
  else {
    // move isn't allowed, remove it if it is in there
    if (SD->MLLr[pos] != -1) {
      if (state == 1) 
	FAremoveFromMLL(SD, 'u', pos);
      else
	FAremoveFromMLL(SD, 'd', pos);
    }
  }
}
inline void EddEast_updateLatPos2(struct SimData *SD, int pos, struct LList *llist) {
  if (LlistLookup(llist, pos))
    // We are already in the list, so we were already processed--
    // don't process it again.
    return;
  LlistAdd(llist, pos);
  EddEast_updateLatPos(SD, pos);
}



int EddEast_consistencyCheck(struct SimData *SD) {
  // first check that all things in MLLr map to the right thing in MLL
  int MLLlocation;
  int moveIndex;    // this is really "pos", but to keep code the same as 
                    // the other event-driven codes i'll use moveIndex.
  int retval=0;
  for (moveIndex=0 ; moveIndex<SD->lattSize ; moveIndex++) {
    int state = SD->lattsite[moveIndex] != S12_EMPTYSITE;
    if (SD->MLLr[moveIndex] != -1) {
      // if it exists in MLLr, it should be at that point in MLL
      if (state == 1 && moveIndex != SD->MLL[SD->MLLr[moveIndex]] ) {
	retval += 1;
	printf("error orjlhc: mi:%d\n", moveIndex);
      } else if (state==0 && moveIndex != SD->MLL_down[SD->MLLr[moveIndex]]) {
	retval += 1;
	printf("error cakrt\n");
      }
    }
  }
  // Now check that all things in the MLL map to the right thing in
  // the MLLr
  for (MLLlocation=0 ; MLLlocation<(SD->lattSize) ; MLLlocation++) {
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
  // Now do the same thing for the MLL_down list:
  // Now check that all things in the MLL map to the right thing in
  // the MLLr
  for (MLLlocation=0 ; MLLlocation<(SD->lattSize) ; MLLlocation++) {
    //if (debugedd)
    //  printf("mlllocation: %d : %d\n", MLLlocation, SD->MLLlen_down);
    if (MLLlocation < SD->MLLlen_down) {
      // if it's less than the list length, then it should be look-up able.
      if (MLLlocation != SD->MLLr[SD->MLL_down[MLLlocation]]) {
	retval += 1;
	printf("error yborkrk\n");
      }
    }
    else {
      // all these greater ones should be blank
      if (SD->MLL_down[MLLlocation] != -1) {
	retval += 1;
	printf("error morkroc\n");
      }
    }
  }
  // Now look and see if everything in MLLr is there if it needs to
  // be... 
  int pos;
  for (pos=0 ; pos<SD->lattSize ; pos++) {
    if (debugedd) printf("--cc: pos: %d\n", pos);
    int atomtype = SD->atomtype[SD->lattsite[pos]];
    if (SD->lattsite[pos] == S12_EMPTYSITE) {
      if (EastNNeighbors(SD, pos) >= SD->inserttype) {
	// It can flip up -- be sure it's in the right list.
        if (SD->MLLr[pos] == -1) {
	  retval += 1;
	  printf("error aroork: not mobile but in MLLr: pos:%d\n", pos);
        }
      } else {
	// be sure that it is not in any of the lookups.
        if (SD->MLLr[pos] != -1) {
	  retval += 1;
	  printf("error vorcktno: not mobile but in MLLr: pos:%d\n", pos);
        }
      }
    } else {
      // So we do have a particle here.  Be sure that it is correct...
      // basically reproduce the logic of the update function.
      
      int isAllowedMove;
      if (EastNNeighbors(SD, pos) >= atomtype)
	isAllowedMove = 1;
      else
	isAllowedMove = 0;
      
      if (isAllowedMove) {
	// We can flip down
	if (SD->MLLr[pos] == -1) {
	  retval += 1;
	  printf("error ycorkon pos:%d\n", pos);
	}
      }
      else {
	// not allowed to flip down
	if (SD->MLLr[pos] != -1) {
	  retval += 1;
	  printf("error xkrgcaon\n");
	}
      }
    }
  }
  return(retval);
}



int EddEast_cycle(struct SimData *SD, double n) {
  if (SD->MLLlen == 0 && SD->MLLlen_down == 0) {
    printf("EddEast_cycle: error, move list length is zero\n");
    exit(12);
  }
  int naccept = 0;
  //struct LList llist; llist.n = 0;

  double maxTime = (double) n;
  double time = SD->MLLextraTime;
  //printf("eddCycle: time:%f maxTime%f\n", time, maxTime);
  int connMax = SD->connMax;

  double c = 1./(1. + exp(SD->beta));
  //printf("beta:%f c:%f\n", SD->beta, c);
  double wUp = 1-(c);
  double wDown = c;

  double wTot = SD->MLLlen_down * wDown   +   SD->MLLlen * wUp;
  if (time == -1.) {
    // pre-move, advance time until the first event.  Otherwise we
    // always end up moving right at time == 0
    double timestep = 1./wTot;
    timestep *= -log(genrand_real3());
    time = timestep;
  }

  while (time < maxTime) {
    int pos;
    double rand = genrand_real2() * wTot;

    if (rand < (SD->MLLlen_down)*wDown) {
      // pick some down spin to flip up
      int i = (int) (rand / wDown);
      pos = SD->MLL_down[i];
      if (debugedd)
	printf("move: flipping down->up at pos:%d\n", pos);
      FAremoveFromMLL(SD, 'd', pos);
      addParticle(SD, pos, SD->inserttype);
      FAaddToMLL(SD, 'u', pos);
      SD->persist[pos] = 1;
    }
    else {
      // pick some up spin to flip down
      int i = (int) ((rand - (SD->MLLlen_down*wDown))  / wUp);
      pos = SD->MLL[i];
      if (debugedd)
	printf("move: flipping up->down at pos:%d\n", pos);
      FAremoveFromMLL(SD, 'u', pos);
      delParticle(SD, pos);
      FAaddToMLL(SD, 'd', pos);
      SD->persist[pos] = 1;
    }

    //llist.n = 0;

    int conni;
    for (conni=0 ; conni < SD->connN[pos] ; conni++) {
      int adjpos = SD->conn[pos*connMax + conni];
      EddEast_updateLatPos(SD, adjpos);
    }
    
    naccept += 1;
    
    // Advance time
    wTot = SD->MLLlen_down * wDown   +   SD->MLLlen * wUp;
    double timestep = 1./wTot;
    timestep *= -log(genrand_real3());  // exponential distribution of times.
                                        // genrand_real3()  -> (0, 1)
    time += timestep;
    //printf("interval: %f\n", (SD->N * SD->connMax) / (double)SD->MLLlen);
    //if (naccept == 1)
    //return(naccept);
  }
  // MLLextraTime is a positive, the amount to increment before our next stop.
  SD->MLLextraTime = time - maxTime;
  //printf("eddCycle: final time:%f maxTime:%f extraTime:%f\n", 
  // time, maxTime, SD->MLLextraTime);

  return(naccept);
}
