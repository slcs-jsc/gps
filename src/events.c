#include "libgps.h"

int main(
  int argc,
  char *argv[]) {

  gps_t *gps;

  FILE *in, *out;

  static double ptmin, ptmax, se[NZ], sz[NZ], w, wmax, wsum2 = 1.0, var;

  static int iarg, ids, idx, iz, sn;

  /* Allocate... */
  ALLOC(gps, gps_t, 1);

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <events.tab> <sens.tab> "
	   "<gps1.nc> [<gps2.nc> ...]");

  /* Create file... */
  printf("Write event data: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = longitude [deg]\n"
	  "# $3 = latitude [deg]\n"
	  "# $4 = minimum perturbation [K]\n"
	  "# $5 = maximum perturbation [K]\n"
	  "# $6 = temperature variance  [K^2]\n\n");

  /* Read vertical sensitivity function... */
  if (argv[3][0] != '-') {
    read_shape(argv[3], sz, se, &sn);
    if (sn > NZ)
      ERRMSG("Too many data points!");
  }

  /* Loop over data files... */
  for (iarg = 4; iarg < argc; iarg++) {

    /* Read gps data... */
    if (!(in = fopen(argv[iarg], "r")))
      continue;
    else {
      fclose(in);
      read_gps(argv[iarg], gps);
    }

    /* Loop over profiles... */
    for (ids = 0; ids < gps->nds; ids++) {

      /* Check tropopause height... */
      if (!gsl_finite(gps->th[ids]))
	continue;

      /* Multiply with vertical sensitivity function... */
      if (argv[3][0] != '-') {
	wmax = wsum2 = 0;
	for (iz = 0; iz < gps->nz[ids]; iz++) {
	  if (gps->z[ids][iz] < sz[0] || gps->z[ids][iz] > sz[sn - 1])
	    w = 0;
	  else {
	    idx = locate_irr(sz, sn, gps->z[ids][iz]);
	    w =
	      LIN(sz[idx], se[idx], sz[idx + 1], se[idx + 1],
		  gps->z[ids][iz]);
	  }
	  if (gsl_finite(gps->t[ids][iz]) && gps->pt[ids][iz]) {
	    gps->pt[ids][iz] *= w;
	    wmax = GSL_MAX(w, wmax);
	    wsum2 += gsl_pow_2(w);
	  }
	}
	for (iz = 0; iz < gps->nz[ids]; iz++)
	  gps->pt[ids][iz] /= wmax;
	wsum2 /= gsl_pow_2(wmax);
      }

      /* Get minimum and maximum perturbation... */
      ptmin = ptmax = var = 0;
      for (iz = 0; iz < gps->nz[ids]; iz++)
	if (gsl_finite(gps->pt[ids][iz])) {
	  ptmin = GSL_MIN(ptmin, gps->pt[ids][iz]);
	  ptmax = GSL_MAX(ptmax, gps->pt[ids][iz]);
	  var += gsl_pow_2(gps->pt[ids][iz]) / wsum2;
	}

      /* Write output... */
      fprintf(out, "%.2f %g %g %g %g %g\n", gps->time[ids],
	      gps->lon[ids][gps->nz[ids] / 2],
	      gps->lat[ids][gps->nz[ids] / 2], ptmin, ptmax, var);
    }
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(gps);

  return EXIT_SUCCESS;
}
