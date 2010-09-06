

inline int energySPM_onSquare(struct SimData *SD,
			      int posA, int posB, int posC, int posD) {
  int ups = (SD->lattsite[posA]!=S12_EMPTYSITE) +
            (SD->lattsite[posB]!=S12_EMPTYSITE) +
            (SD->lattsite[posC]!=S12_EMPTYSITE) +
            (SD->lattsite[posD]!=S12_EMPTYSITE);
  // Use the formula:  (-1/2) * (s_1*s_2*s_3*s_4 - 1)
  // So there is an implicit factor of 1/2
  // Even number of spins makes energy 0
  // Odd number of spins makes energy -1
  if ((ups % 2) == 0)
    return 0;
  else
    return -1;
}

inline double energySPM_posLocal(struct SimData *SD, int pos) {


  int N  = SD->conn[pos*4+0];
  int NE = SD->conn[N  *4+1];
  int E  = SD->conn[pos*4+1];

  if (errorcheck) {
    if (NE != SD->conn[E*4+0]) {
      printf("Error tnsako6451ntko\n");
      exit(75);
    }
  }

  return SD->hardness * energySPM_onSquare(SD, pos, N, NE, E);
}

inline double energySPM_pos(struct SimData *SD, int pos) {
  int N  = SD->conn[pos*4+0];
  int NE = SD->conn[N  *4+1];
  int E  = SD->conn[pos*4+1];
  int ES = SD->conn[E  *4+2];
  int S  = SD->conn[pos*4+2];
  int SW = SD->conn[S  *4+3];
  int W  = SD->conn[pos*4+3];
  int WN = SD->conn[W  *4+0];

  return SD->hardness * (energySPM_onSquare(SD, pos, N, NE, E) +
			 energySPM_onSquare(SD, pos, E, ES, S) +
			 energySPM_onSquare(SD, pos, S, SW, W) +
			 energySPM_onSquare(SD, pos, W, WN, N));
}




