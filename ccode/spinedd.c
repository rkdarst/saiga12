/* Richard Darst, September 2010 */

/* These functions were coppied from fredricksonandersen.c since they
 *  are used multiple places */


inline void FAaddToMLL(struct SimData *SD, char which, int pos) {
  int *MLL, *MLLlen;
  if (errorcheck)
    if (which != 'u' && which != 'd')
      printf("error: which is not 'u' or 'd': %c\n", which);
  if (which == 'u') {
    if(debugedd) 
      printf("adding to up list: %c %d\n", which, pos);
    MLL    =   SD->MLL;
    //MLLr   =   SD->MLLr;
    MLLlen = &(SD->MLLlen);
  }
  else {
    if(debugedd) 
      printf("adding to down list: %c %d\n", which, pos);
    MLL    =   SD->MLL_down;
    //MLLr   =   SD->MLLr_down;
    MLLlen = &(SD->MLLlen_down);
  }
  if (errorcheck && SD->MLLr[pos] != -1) {
    printf("error rhoeacurk %c %d \n", which, pos);
  }
  MLL[*MLLlen] = pos;
  SD->MLLr[pos] = *MLLlen;
  *MLLlen += 1;
}
inline void FAremoveFromMLL(struct SimData *SD, char which, int pos) {
  int *MLL, *MLLlen;
  if (errorcheck)
    if (which != 'u' && which != 'd')
      printf("error: which is not 'u' or 'd': %c\n", which);
    
  if (which == 'u') {
    if (errorcheck)
      if (SD->MLLr[pos] == -1)
	printf("error: trying to remove from up list but is not up:%d\n", pos);
    MLL    =   SD->MLL;
    //MLLr   =   SD->MLLr;
    MLLlen = &(SD->MLLlen);
  }
  else {
    if (errorcheck)
      if (SD->MLLr[pos] == -1)
	printf("error: trying to remove from down list but is not dn:%d\n", pos);
    MLL    =   SD->MLL_down;
    //MLLr   =   SD->MLLr_down;
    MLLlen = &(SD->MLLlen_down);
  }
  int mllLocation = SD->MLLr[pos];
  if (mllLocation == (*MLLlen)-1 ) {
    MLL[mllLocation] = -1;
    SD->MLLr[pos] = -1;
  } else {
    MLL[mllLocation] = MLL[(*MLLlen)-1];
    SD->MLLr[MLL[(*MLLlen)-1]] = mllLocation;

    SD->MLLr[pos] = -1;
    MLL[(*MLLlen)-1] = -1;
  }
  *MLLlen -= 1 ;
}
