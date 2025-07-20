/*
  This file is part of the GPS Code Collection.
  
  the GPS Code Collections is free software: you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.
  
  The GPS Code Collection is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with the GPS Code Collection. If not, see
  <http://www.gnu.org/licenses/>.
  
  Copyright (C) 2019-2025 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Extract gravity wave perturbations.
*/

#include "libgps.h"

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  gps_t *gps;

  char metbase[LEN];

  /* Allocate... */
  ALLOC(gps, gps_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <out.nc> <gps1.nc> [<gps2.nc> ...]");

  /* Get control parameters... */
  const double dt_met = scan_ctl(argc, argv, "DT_MET", -1, "21600", NULL);
  const double gauss_dx = scan_ctl(argc, argv, "GAUSS_DX", -1, "-999", NULL);
  const double gauss_dy = scan_ctl(argc, argv, "GAUSS_DY", -1, "-999", NULL);
  const double grid_zmin = scan_ctl(argc, argv, "GRID_ZMIN", -1, "0", NULL);
  const double grid_zmax = scan_ctl(argc, argv, "GRID_ZMAX", -1, "40", NULL);
  const int grid_nz = (int) scan_ctl(argc, argv, "GRID_NZ", -1, "-1", NULL);
  const double ham_dz = scan_ctl(argc, argv, "HAM_DZ", -1, "-999", NULL);
  const double ham_dz2 = scan_ctl(argc, argv, "HAM_DZ2", -1, "-999", NULL);
  scan_ctl(argc, argv, "METBASE", -1, "", metbase);
  const int poly_dim = (int) scan_ctl(argc, argv, "POLY_DIM", -1, "5", NULL);
  const double poly_zmin = scan_ctl(argc, argv, "POLY_ZMIN", -1, "0", NULL);
  const double poly_zmax = scan_ctl(argc, argv, "POLY_ZMAX", -1, "40", NULL);

  /* Read individual GPS-RO data files... */
  for (int iarg = 3; iarg < argc; iarg++) {
    FILE *in;
    if (!(in = fopen(argv[iarg], "r")))
      continue;
    else {
      fclose(in);
      read_gps_prof(argv[iarg], gps);
    }
  }

  /* Check number of profiles... */
  if (gps->nds <= 0)
    ERRMSG("No profiles found!");

  /* Grid profile... */
  if (grid_nz > 0)
    grid_gps(gps, grid_zmin, grid_zmax, grid_nz);

  /* Get tropopause... */
  tropopause(gps);

  /* Get perturbations from horizontal Gaussian mean... */
  if (gauss_dx > 0 && gauss_dy > 0)
    gauss(gps, gauss_dx, gauss_dy);

  /* Get perturbations from vertical Hamming filter... */
  if (ham_dz > 0)
    hamming_low_pass(gps, ham_dz);

  /* Use vertical Hamming filter to reduce noise... */
  if (ham_dz2 > 0)
    hamming_high_pass(gps, ham_dz2);

  /* Use meteo data for detrending... */
  if (metbase[0] != '-')
    detrend_met(gps, metbase, dt_met);

  /* Remove polynomial fit from perturbation profile... */
  if (poly_dim > 0)
    poly(gps, poly_dim, poly_zmin, poly_zmax);

  /* Write GPS-RO data file... */
  write_gps(argv[2], gps);

  /* Free... */
  free(gps);

  return EXIT_SUCCESS;
}
