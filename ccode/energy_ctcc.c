/* Richard Darst, May 2009 */


inline double energyCTCC_posLocal(struct SimData *SD, int pos) {
  if (SD->lattsite[pos] == S12_EMPTYSITE)
    return 0;
  //int type = atomType(SD, pos);
  //int excessneighbors = SD->nneighbors[pos] - type;
  //printf("excessneighbors: %d\n", excessneighbors);
  int orient = SD->orient[pos];
  int neighPos = SD->conn[SD->connMax*pos + orient];
  if (SD->lattsite[neighPos] == S12_EMPTYSITE)
    return 0.;
  else
    return SD->hardness;
}

inline double energyCTCC_pos(struct SimData *SD, int pos) {
  if (SD->lattsite[pos] == S12_EMPTYSITE)
    return 0;
  int i_conn;
  double E = energyCTCC_posLocal(SD, pos);
  for (i_conn=0; i_conn<SD->connN[pos] ; i_conn++) {
    int adjPos = SD->conn[pos*SD->connMax + i_conn];
    int reverseconn = (i_conn+(SD->connMax/2)) % SD->connMax;
    if(SD->orient[adjPos] == reverseconn) {
      E += SD->hardness;
    }
  }
  return E;
}
