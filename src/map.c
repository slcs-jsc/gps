#include "libgps.h"

int main(
  int argc,
  char *argv[]) {

  gps_t *gps;

  FILE *out;

  double z;

  int ids, iz;

  /* Allocate... */
  ALLOC(gps, gps_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <gps.nc> <map.tab>");

  /* Get control parameters... */
  z = scan_ctl(argc, argv, "Z", -1, "20", NULL);

  /* Read gps data... */
  read_gps(argv[2], gps);

  /* Create output file... */
  printf("Write map data: %s\n", argv[3]);
  if (!(out = fopen(argv[3], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [sec]\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = longitude [deg]\n"
	  "# $4 = latitude [deg]\n"
	  "# $5 = pressure [hPa]\n"
	  "# $6 = temperature [K]\n"
	  "# $7 = water vapor vmr [ppm]\n"
	  "# $8 = temperature perturbation [K]\n"
	  "# $9 = tropopause height [km]\n\n");

  /* Write data... */
  for (ids = 0; ids < gps->nds; ids++)
    for (iz = 0; iz < gps->nz[ids]; iz++)
      if (fabs(gps->z[ids][iz] - z) < 0.01) {
	fprintf(out, "%.2f %g %g %g %g %g %g %g %g\n",
		gps->time[ids], gps->z[ids][iz], gps->lon[ids][iz],
		gps->lat[ids][iz], gps->p[ids][iz], gps->t[ids][iz],
		gps->wv[ids][iz], gps->pt[ids][iz], gps->th[ids]);
	break;
      }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(gps);

  return EXIT_SUCCESS;
}
