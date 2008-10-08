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
