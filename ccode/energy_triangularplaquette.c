

inline double energyTPM_onSquare(struct SimData *SD,
				 int posA, int posB, int posC) {
  int ups = (SD->lattsite[posA]!=S12_EMPTYSITE) +
            (SD->lattsite[posB]!=S12_EMPTYSITE) +
            (SD->lattsite[posC]!=S12_EMPTYSITE);
  // Use the formula:  (-1/2) * (s_1*s_2*s_3)
  // So there is an implicit factor of 1/2
  // Even number of spins makes energy -.5
  // Odd number of spins makes energy +.5
  if ((ups % 2) == 0)
    return +.5;  // should be +.5...
  else
    return -.5;
}

inline double energyTPM_posLocal(struct SimData *SD, int pos) {


  int N  = SD->conn[pos*4+0];
  //int NE = SD->conn[N  *4+1];
  int E  = SD->conn[pos*4+1];

  /* if (errorcheck) { */
  /*   if (NE != SD->conn[E*4+0]) { */
  /*     printf("Error tnsako6451ntko\n"); */
  /*     exit(75); */
  /*   } */
  /* } */

  return SD->hardness * energyTPM_onSquare(SD, pos, N, E);
}

inline double energyTPM_pos(struct SimData *SD, int pos) {
  int N  = SD->conn[pos*4+0];
  //int NE = SD->conn[N  *4+1];
  int E  = SD->conn[pos*4+1];
  int ES = SD->conn[E  *4+2];
  int S  = SD->conn[pos*4+2];
  //int SW = SD->conn[S  *4+3];
  int W  = SD->conn[pos*4+3];
  int WN = SD->conn[W  *4+0];

  return SD->hardness * (energyTPM_onSquare(SD, pos, N,  E) +
			 energyTPM_onSquare(SD, pos, ES, S) +
			 //energySPM_onSquare(SD, S, SW, W) +
			 energyTPM_onSquare(SD, pos, W, WN));
}




