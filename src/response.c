#include "libgps.h"

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  gps_t *gps;

  FILE *out;

  double lz, ptmax[NZ], var[NZ], w, wmax, se[NZ], sz[NZ], t0 = 10.0,
    grid_zmin, grid_zmax, ham_dz, ham_dz2;

  int idx, iphi, iz, iz2, nphi = 360, sn = 0, grid_nz;

  /* Allocate... */
  ALLOC(gps, gps_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <sens.tab> <response.tab>");

  /* Get control parameters... */
  grid_zmin = scan_ctl(argc, argv, "GRID_ZMIN", -1, "0", NULL);
  grid_zmax = scan_ctl(argc, argv, "GRID_ZMAX", -1, "60", NULL);
  grid_nz = (int) scan_ctl(argc, argv, "GRID_NZ", -1, "601", NULL);
  ham_dz = scan_ctl(argc, argv, "HAM_DZ", -1, "6.0", NULL);
  ham_dz2 = scan_ctl(argc, argv, "HAM_DZ2", -1, "0.4", NULL);

  /* Read vertical sensitivity function... */
  if (argv[2][0] != '-') {
    read_shape(argv[2], sz, se, &sn);
    if (sn > NZ)
      ERRMSG("Too many data points!");
  }

  /* Create output file... */
  printf("Write response data: %s\n", argv[3]);
  if (!(out = fopen(argv[3], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = vertical wavelength [km]\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = response (amplitude) [%%]\n"
	  "# $4 = response (variance) [%%]\n");

  /* Create profile... */
  gps->nds = 1;
  gps->nz[0] = grid_nz;
  for (iz = 0; iz < gps->nz[0]; iz++)
    gps->z[0][iz] =
      LIN(0.0, grid_zmin, grid_nz - 1.0, grid_zmax, (double) iz);

  /* Loop over vertical wavelength... */
  for (lz = 0.1; lz <= 20.0; lz += 0.1) {

    /* Write info... */
    printf("Calculate %g km...\n", lz);

    /* Initialize... */
    for (iz = 0; iz < gps->nz[0]; iz++)
      ptmax[iz] = var[iz] = 0;

    /* Loop over phase... */
    for (iphi = 0; iphi < nphi; iphi++) {

      /* Create profile... */
      for (iz = 0; iz < gps->nz[0]; iz++)
	gps->t[0][iz] = 250.0 + t0 * sin(2. * M_PI / lz * gps->z[0][iz]
					 +
					 2. * M_PI * (double) iphi /
					 (double) nphi);

      /* Get perturbations from vertical Hamming filter... */
      if (ham_dz > 0)
	hamming_low_pass(gps, ham_dz);

      /* Use vertical Hamming filter to reduce noise... */
      if (ham_dz2 > 0)
	hamming_high_pass(gps, ham_dz2);

      /* Multiply with vertical sensitivity function... */
      if (argv[2][0] != '-') {
	wmax = 0;
	for (iz2 = 0; iz2 < gps->nz[0]; iz2++) {
	  if (gps->z[0][iz2] < sz[0] || gps->z[0][iz2] > sz[sn - 1])
	    w = 0;
	  else {
	    idx = locate_irr(sz, sn, gps->z[0][iz2]);
	    w =
	      LIN(sz[idx], se[idx], sz[idx + 1], se[idx + 1], gps->z[0][iz2]);
	  }
	  gps->pt[0][iz2] *= w;
	  wmax = GSL_MAX(w, wmax);
	}
	if (wmax > 0)
	  for (iz2 = 0; iz2 < gps->nz[0]; iz2++)
	    gps->pt[0][iz2] /= wmax;
      }

      /* Get response... */
      for (iz = 0; iz < gps->nz[0]; iz++) {
	ptmax[iz] = GSL_MAX(ptmax[iz], gps->pt[0][iz]);
	var[iz] += gsl_pow_2(gps->pt[0][iz]) / nphi;
      }
    }

    /* Write output... */
    fprintf(out, "\n");
    for (iz = 0; iz < gps->nz[0]; iz++)
      fprintf(out, "%g %g %g %g\n", lz, gps->z[0][iz], ptmax[iz] / t0,
	      2.0 * var[iz] / gsl_pow_2(t0));
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(gps);

  return EXIT_SUCCESS;
}
