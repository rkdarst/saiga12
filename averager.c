
struct averager_dat {
  long n;
  double mean;
  double M2
};

void averager_add (struct averager_dat *dat, double val) {
  long n = dat->n;
  double delta = val - dat->mean;
  float mean = dat->mean + delta/n;
  double M2 = dat->M2 + delta(val - mean);
  dat->n = n;
  dat->mean = mean;
  dat->M2 = M2;
}
double averager_mean (struct averager_dat *dat) {
  return(dat->mean);
}
double averager_std (struct averager_dat *dat) {
  return(sqrt(dat->M2 / dat->n));
}
double averager_var (struct averager_dat *dat) {
  return(dat->M2 / dat->n);
}
double averager_stdsample (struct averager_dat *dat) {
  return(sqrt(dat->M2 / (dat->n-1)));
}
double averager_varsample (struct averager_dat *dat) {
  return(dat->M2 / (dat->n-1));
}

