/* Richard Darst, October 2008 */

inline double energyBM_posLocal(struct SimData *SD, int pos) {
  if (SD->lattsite[pos] == S12_EMPTYSITE)
    return 0;
  int type = atomType(SD, pos);
  int excessneighbors = SD->nneighbors[pos] - type;
  //printf("excessneighbors: %d\n", excessneighbors);
  if (excessneighbors <= 0) {
    return 0.;
  }
  else {
    return excessneighbors * SD->hardness;
  }
}

inline double energyBM_pos(struct SimData *SD, int pos) {
  int i_conn;
  double E = energyBM_posLocal(SD, pos);
  for (i_conn=0; i_conn<SD->connN[pos] ; i_conn++) {
    int adjPos = SD->conn[pos*SD->connMax + i_conn];
    E += energyBM_posLocal(SD, adjPos);
  }
  return E;
}


/*
 *  Configs are heavily penalized for having zero neighbors.
 */
inline double energyBMnotzero_posLocal(struct SimData *SD, int pos) {
  if (SD->lattsite[pos] == S12_EMPTYSITE)
    return 0;
  int type = atomType(SD, pos);
  int excessneighbors = SD->nneighbors[pos] - type;
  //printf("excessneighbors: %d\n", excessneighbors);
  if (SD->nneighbors[pos] == 0)
    return (SD->hardness * 10.) ;
  if (excessneighbors <= 0)
    return 0.;
  else
    return excessneighbors * SD->hardness;
}

inline double energyBMnotzero_pos(struct SimData *SD, int pos) {
  int i_conn;
  double E = energyBMnotzero_posLocal(SD, pos);
  for (i_conn=0; i_conn<SD->connN[pos] ; i_conn++) {
    int adjPos = SD->conn[pos*SD->connMax + i_conn];
    E += energyBMnotzero_posLocal(SD, adjPos);
  }
  return E;
}


/*
 *  Packing violations are much more penalized much more for type 1
 *  particles.
 */
inline double energyBMimmobile1_posLocal(struct SimData *SD, int pos) {
  if (SD->lattsite[pos] == S12_EMPTYSITE)
    return 0;
  int type = atomType(SD, pos);
  int excessneighbors = SD->nneighbors[pos] - type;
  //printf("excessneighbors: %d\n", excessneighbors);
  if (excessneighbors <= 0)
    return 0.;
  else if (type == 1)
    return (excessneighbors * SD->hardness * 10);
  else
    return excessneighbors * SD->hardness;

}

inline double energyBMimmobile1_pos(struct SimData *SD, int pos) {
  int i_conn;
  double E = energyBMimmobile1_posLocal(SD, pos);
  for (i_conn=0; i_conn<SD->connN[pos] ; i_conn++) {
    int adjPos = SD->conn[pos*SD->connMax + i_conn];
    E += energyBMimmobile1_posLocal(SD, adjPos);
  }
  return E;
}

/*
 *  Packing violations are penalized infinitely for type 1 particles,
 *  using only 'hardness' for other types.
 */
inline double energyBMimmobile1b_posLocal(struct SimData *SD, int pos) {
  if (SD->lattsite[pos] == S12_EMPTYSITE)
    return 0;
  int type = atomType(SD, pos);
  int excessneighbors = SD->nneighbors[pos] - type;
  //printf("excessneighbors: %d\n", excessneighbors);
  if (excessneighbors <= 0)
    return 0.;
  else if (type == 1)
    return 1./0.;
  else
    return excessneighbors * SD->hardness;
}

inline double energyBMimmobile1b_pos(struct SimData *SD, int pos) {
  int i_conn;
  double E = energyBMimmobile1_posLocal(SD, pos);
  for (i_conn=0; i_conn<SD->connN[pos] ; i_conn++) {
    int adjPos = SD->conn[pos*SD->connMax + i_conn];
    E += energyBMimmobile1_posLocal(SD, adjPos);
  }
  return E;
}
