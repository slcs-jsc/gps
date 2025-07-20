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
  Estimate gravity wave variances.
*/

#include "libgps.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum number of longitudes. */
#define GX 360

/* Maximum number of latitudes. */
#define GY 180

/* Maximum number of altitudes. */
#define GZ 50

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  gps_t *gps;

  FILE *out;

  static double lon0, lon1, lat0, lat1, z0, z1, mean[GX][GY][GZ],
    min[GX][GY][GZ], max[GX][GY][GZ], var[GX][GY][GZ],
    mtime[GX][GY], glon[GX], glat[GY], gz[GZ], thmean[GX][GY],
    tmean[GX][GY][GZ], twmean[GX][GY], se[NZ], sz[NZ], tw, w, wmax, wsum;

  static int iarg, ids, idx, ix, iy, iz, iz2,
    nx, ny, nz, n[GX][GY][GZ], np[GX][GY], sn;

  /* Allocate... */
  ALLOC(gps, gps_t, 1);

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <var.tab> <sens.tab> "
	   "<gps1.nc> [<gps2.nx> ...]");

  /* Get control parameters... */
  z0 = scan_ctl(argc, argv, "Z0", -1, "0", NULL);
  z1 = scan_ctl(argc, argv, "Z1", -1, "60", NULL);
  nz = (int) scan_ctl(argc, argv, "NZ", -1, "6", NULL);
  lon0 = scan_ctl(argc, argv, "LON0", -1, "-180", NULL);
  lon1 = scan_ctl(argc, argv, "LON1", -1, "180", NULL);
  nx = (int) scan_ctl(argc, argv, "NX", -1, "72", NULL);
  lat0 = scan_ctl(argc, argv, "LAT0", -1, "-90", NULL);
  lat1 = scan_ctl(argc, argv, "LAT1", -1, "90", NULL);
  ny = (int) scan_ctl(argc, argv, "NY", -1, "36", NULL);

  /* Check grid dimensions... */
  if (nx < 1 || nx > GX)
    ERRMSG("Set 1 <= GX <= MAX!");
  if (ny < 1 || ny > GY)
    ERRMSG("Set 1 <= GY <= MAX!");
  if (nz < 1 || nz > GZ)
    ERRMSG("Set 1 <= GZ <= MAX!");

  /* Read vertical sensitivity function... */
  if (argv[3][0] != '-') {
    read_shape(argv[3], sz, se, &sn);
    if (sn > NZ)
      ERRMSG("Too many data points!");
  }

  /* Loop over data files... */
  for (iarg = 4; iarg < argc; iarg++) {

    /* Read gps data... */
    FILE *in;
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
	tw = wsum = wmax = 0;
	for (iz2 = 0; iz2 < gps->nz[ids]; iz2++) {
	  if (gps->z[ids][iz2] < sz[0] || gps->z[ids][iz2] > sz[sn - 1])
	    w = 0;
	  else {
	    idx = locate_irr(sz, sn, gps->z[ids][iz2]);
	    w =
	      LIN(sz[idx], se[idx], sz[idx + 1], se[idx + 1],
		  gps->z[ids][iz2]);
	  }
	  if (gsl_finite(gps->t[ids][iz2]) && gps->pt[ids][iz2]) {
	    tw += w * gps->t[ids][iz2];
	    wsum += w;
	    gps->pt[ids][iz2] *= w;
	    wmax = GSL_MAX(w, wmax);
	  }
	}
	tw /= wsum;
	for (iz2 = 0; iz2 < gps->nz[ids]; iz2++)
	  gps->pt[ids][iz2] /= wmax;
      }

      /* Get grid indices... */
      ix = (int) ((gps->lon[ids][gps->nz[ids] / 2] - lon0)
		  / (lon1 - lon0) * (double) nx);
      iy = (int) ((gps->lat[ids][gps->nz[ids] / 2] - lat0)
		  / (lat1 - lat0) * (double) ny);
      if (ix < 0 || ix >= nx || iy < 0 || iy >= ny)
	continue;

      /* Get mean time and tropopause height... */
      mtime[ix][iy] += gps->time[ids];
      thmean[ix][iy] += gps->th[ids];
      twmean[ix][iy] += tw;
      np[ix][iy]++;

      /* Loop over altitudes... */
      for (iz2 = 0; iz2 < gps->nz[ids]; iz2++) {

	/* Get grid indices... */
	iz = (int) ((gps->z[ids][iz2] - z0)
		    / (z1 - z0) * (double) nz);
	if (iz < 0 || iz >= nz)
	  continue;

	/* Check data... */
	if (!gsl_finite(gps->t[ids][iz2]) || !gsl_finite(gps->pt[ids][iz2]))
	  continue;

	/* Get statistics of perturbations... */
	tmean[ix][iy][iz] += gps->t[ids][iz2];
	mean[ix][iy][iz] += gps->pt[ids][iz2];
	var[ix][iy][iz] += gsl_pow_2(gps->pt[ids][iz2]);
	max[ix][iy][iz] = GSL_MAX(max[ix][iy][iz], gps->pt[ids][iz2]);
	min[ix][iy][iz] = GSL_MIN(min[ix][iy][iz], gps->pt[ids][iz2]);
	n[ix][iy][iz]++;
      }
    }
  }

  /* Analyze results... */
  for (ix = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++) {

      /* Get mean time and tropopause height... */
      if (np[ix][iy] > 0) {
	mtime[ix][iy] /= (double) np[ix][iy];
	thmean[ix][iy] /= (double) np[ix][iy];
	twmean[ix][iy] /= (double) np[ix][iy];
      } else {
	mtime[ix][iy] = GSL_NAN;
	thmean[ix][iy] = GSL_NAN;
	twmean[ix][iy] = GSL_NAN;
      }

      /* Loop over altitudes... */
      for (iz = 0; iz < nz; iz++) {

	/* Get geolocation... */
	gz[iz] = z0 + (iz + 0.5) / (double) nz *(
  z1 - z0);
	glon[ix] = lon0 + (ix + 0.5) / (double) nx *(
  lon1 - lon0);
	glat[iy] = lat0 + (iy + 0.5) / (double) ny *(
  lat1 - lat0);

	/* Get mean perturbation and variance... */
	if (n[ix][iy][iz] > 0) {
	  tmean[ix][iy][iz]
	    /= (double) n[ix][iy][iz];
	  mean[ix][iy][iz]
	    /= (double) n[ix][iy][iz];
	  var[ix][iy][iz]
	    = var[ix][iy][iz] / (double) n[ix][iy][iz]
	    - gsl_pow_2(mean[ix][iy][iz]);
	} else {
	  tmean[ix][iy][iz] = GSL_NAN;
	  mean[ix][iy][iz] = GSL_NAN;
	  var[ix][iy][iz] = GSL_NAN;
	  min[ix][iy][iz] = GSL_NAN;
	  max[ix][iy][iz] = GSL_NAN;
	}
      }
    }

  /* Create file... */
  printf("Write variance statistics: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = longitude [deg]\n"
	  "# $4 = latitude [deg]\n"
	  "# $5 = number of profiles\n"
	  "# $6 = number of data points\n"
	  "# $7 = mean perturbation [K]\n"
	  "# $8 = minimum perturbation [K]\n"
	  "# $9 = maximum perturbation [K]\n"
	  "# $10 = variance [K^2]\n"
	  "# $11 = mean temperature [K]\n"
	  "# $12 = mean weighted temperature [K]\n"
	  "# $13 = mean tropopause height [km]\n");

  /* Write results... */
  for (iz = 0; iz < nz; iz++) {
    for (iy = 0; iy < ny; iy++) {
      if (iy == 0 || nx > 1)
	fprintf(out, "\n");
      for (ix = 0; ix < nx; ix++)
	fprintf(out, "%.2f %g %g %g %d %d %g %g %g %g %g %g %g\n",
		mtime[ix][iy], gz[iz], glon[ix], glat[iy],
		np[ix][iy], n[ix][iy][iz], mean[ix][iy][iz],
		min[ix][iy][iz], max[ix][iy][iz], var[ix][iy][iz],
		tmean[ix][iy][iz], twmean[ix][iy], thmean[ix][iy]);
    }
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(gps);

  return EXIT_SUCCESS;
}
