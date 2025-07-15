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
  GPS Code Collection library definitions.
*/

#include "libgps.h"

/*****************************************************************************/

void add_var(
  int ncid,
  const char *varname,
  const char *unit,
  const char *longname,
  int type,
  int dimid[],
  int *varid,
  int ndims) {

  double dp = GSL_NAN;

  /* Define variable... */
  NC(nc_def_var(ncid, varname, type, ndims, dimid, varid));

  /* Set long name... */
  NC(nc_put_att_text(ncid, *varid, "long_name", strlen(longname), longname));

  /* Set units... */
  NC(nc_put_att_text(ncid, *varid, "units", strlen(unit), unit));

  /* Set fill value... */
  NC(nc_put_att_double(ncid, *varid, "_FillValue", type, 1, &dp));
}

/*****************************************************************************/

void detrend_met(
  gps_t *gps,
  char *metbase,
  double dt_met) {

  met_t *met0, *met1;

  double t;

  int ids, iz;

  /* Allocate... */
  ALLOC(met0, met_t, 1);
  ALLOC(met1, met_t, 1);

  /* Loop over profiles... */
  for (ids = 0; ids < gps->nds; ids++) {

    /* Loop over altitudes... */
    for (iz = 0; iz < gps->nz[ids]; iz++) {

      /* Get meteorological data... */
      get_met(metbase, dt_met, gps->time[ids], met0, met1);

      /* Interpolate meteorological data... */
      intpol_met_time(met0, met1, gps->time[ids], gps->p[ids][iz],
		      gps->lon[ids][iz], gps->lat[ids][iz], &t);

      /* Set perturbation... */
      gps->pt[ids][iz] = gps->t[ids][iz] - t;
    }
  }

  /* Free... */
  free(met0);
  free(met1);
}

/*****************************************************************************/

void gauss(
  gps_t *gps,
  double dx,
  double dy) {

  double dlat, dlon, w, wsum;

  int ids, ids2, iz;

  /* Loop over profiles... */
  for (ids = 0; ids < gps->nds; ids++) {

    /* Initialize... */
    wsum = 0;
    for (iz = 0; iz < gps->nz[ids]; iz++)
      gps->pt[ids][iz] = 0;

    /* Calculate lon-lat standard deviations... */
    dlat = dx * 180. / (M_PI * RE) / 2.3548;
    dlon = dy * 180. / 2.3548
      / (M_PI * RE * cos(gps->lat[ids][gps->nz[ids] / 2] * M_PI / 180.));

    /* Calculate mean temperature... */
    for (ids2 = 0; ids2 < gps->nds; ids2++) {
      w = exp(-0.5 * gsl_pow_2((gps->lon[ids][gps->nz[ids] / 2]
				- gps->lon[ids2][gps->nz[ids2] / 2]) / dlon)
	      - 0.5 * gsl_pow_2((gps->lat[ids][gps->nz[ids] / 2]
				 -
				 gps->lat[ids2][gps->nz[ids2] / 2]) / dlat));
      wsum += w;
      for (iz = 0; iz < gps->nz[ids]; iz++)
	gps->pt[ids][iz] += w * gps->t[ids2][iz];
    }

    /* Normalize... */
    if (wsum > 0)
      for (iz = 0; iz < gps->nz[ids]; iz++)
	gps->pt[ids][iz] = gps->t[ids][iz] - gps->pt[ids][iz] / wsum;
  }
}

/*****************************************************************************/

void grid_gps(
  gps_t *gps,
  double zmin,
  double zmax,
  int nz) {

  double lat[NZ], lon[NZ], p[NZ], pt[NZ], t[NZ], wv[NZ], z[NZ];

  int ids, iz, iz2;

  /* Check number of altitudes... */
  if (nz > NZ)
    ERRMSG("Too many altitudes!");

  /* Loop over profiles... */
  for (ids = 0; ids < gps->nds; ids++) {

    /* Loop over altitudes... */
    for (iz = 0; iz < nz; iz++) {

      /* Set altitude... */
      z[iz] = LIN(0.0, zmin, nz - 1.0, zmax, (double) iz);

      /* Get index... */
      iz2 = locate_irr(gps->z[ids], gps->nz[ids], z[iz]);

      /* Interpolate... */
      lon[iz] = LIN(gps->z[ids][iz2], gps->lon[ids][iz2],
		    gps->z[ids][iz2 + 1], gps->lon[ids][iz2 + 1], z[iz]);
      lat[iz] = LIN(gps->z[ids][iz2], gps->lat[ids][iz2],
		    gps->z[ids][iz2 + 1], gps->lat[ids][iz2 + 1], z[iz]);
      p[iz] = LIN(gps->z[ids][iz2], gps->p[ids][iz2],
		  gps->z[ids][iz2 + 1], gps->p[ids][iz2 + 1], z[iz]);
      t[iz] = LIN(gps->z[ids][iz2], gps->t[ids][iz2],
		  gps->z[ids][iz2 + 1], gps->t[ids][iz2 + 1], z[iz]);
      wv[iz] = LIN(gps->z[ids][iz2], gps->wv[ids][iz2],
		   gps->z[ids][iz2 + 1], gps->wv[ids][iz2 + 1], z[iz]);
      pt[iz] = LIN(gps->z[ids][iz2], gps->pt[ids][iz2],
		   gps->z[ids][iz2 + 1], gps->pt[ids][iz2 + 1], z[iz]);
    }

    /* Copy data... */
    gps->nz[ids] = nz;
    for (iz = 0; iz < nz; iz++) {
      gps->z[ids][iz] = z[iz];
      gps->lon[ids][iz] = lon[iz];
      gps->lat[ids][iz] = lat[iz];
      gps->p[ids][iz] = p[iz];
      gps->t[ids][iz] = t[iz];
      gps->wv[ids][iz] = wv[iz];
      gps->pt[ids][iz] = pt[iz];
    }
  }
}

/*****************************************************************************/

void get_met(
  char *metbase,
  double dt_met,
  double t,
  met_t *met0,
  met_t *met1) {

  char filename[LEN];

  static int init;

  /* Init... */
  if (!init) {
    init = 1;

    get_met_help(t, -1, metbase, dt_met, filename);
    read_met(filename, met0);

    get_met_help(t + 1.0, 1, metbase, dt_met, filename);
    read_met(filename, met1);
  }

  /* Read new data... */
  if (t > met1->time) {
    memcpy(met0, met1, sizeof(met_t));
    get_met_help(t, 1, metbase, dt_met, filename);
    read_met(filename, met1);
  }
}

/*****************************************************************************/

void get_met_help(
  double t,
  int direct,
  char *metbase,
  double dt_met,
  char *filename) {

  double t6, r;

  int year, mon, day, hour, min, sec;

  /* Round time to fixed intervals... */
  if (direct == -1)
    t6 = floor(t / dt_met) * dt_met;
  else
    t6 = ceil(t / dt_met) * dt_met;

  /* Decode time... */
  jsec2time(t6, &year, &mon, &day, &hour, &min, &sec, &r);

  /* Set filename... */
  sprintf(filename, "%s_%d_%02d_%02d_%02d.nc", metbase, year, mon, day, hour);
}

/*****************************************************************************/

void intpol_met_3d(
  float array[EX][EY][EP],
  int ip,
  int ix,
  int iy,
  double wp,
  double wx,
  double wy,
  double *var) {

  double aux00, aux01, aux10, aux11;

  /* Interpolate vertically... */
  aux00 = wp * (array[ix][iy][ip] - array[ix][iy][ip + 1])
    + array[ix][iy][ip + 1];
  aux01 = wp * (array[ix][iy + 1][ip] - array[ix][iy + 1][ip + 1])
    + array[ix][iy + 1][ip + 1];
  aux10 = wp * (array[ix + 1][iy][ip] - array[ix + 1][iy][ip + 1])
    + array[ix + 1][iy][ip + 1];
  aux11 = wp * (array[ix + 1][iy + 1][ip] - array[ix + 1][iy + 1][ip + 1])
    + array[ix + 1][iy + 1][ip + 1];

  /* Interpolate horizontally... */
  aux00 = wy * (aux00 - aux01) + aux01;
  aux11 = wy * (aux10 - aux11) + aux11;
  *var = wx * (aux00 - aux11) + aux11;
}

/*****************************************************************************/

void intpol_met_space(
  met_t *met,
  double p,
  double lon,
  double lat,
  double *t) {

  double wp, wx, wy;

  int ip, ix, iy;

  /* Check longitude... */
  if (met->lon[met->nx - 1] > 180 && lon < 0)
    lon += 360;

  /* Get indices... */
  ip = locate_irr(met->p, met->np, p);
  ix = locate_reg(met->lon, met->nx, lon);
  iy = locate_reg(met->lat, met->ny, lat);

  /* Get weights... */
  wp = (met->p[ip + 1] - p) / (met->p[ip + 1] - met->p[ip]);
  wx = (met->lon[ix + 1] - lon) / (met->lon[ix + 1] - met->lon[ix]);
  wy = (met->lat[iy + 1] - lat) / (met->lat[iy + 1] - met->lat[iy]);

  /* Interpolate... */
  intpol_met_3d(met->t, ip, ix, iy, wp, wx, wy, t);
}

/*****************************************************************************/

void intpol_met_time(
  met_t *met0,
  met_t *met1,
  double ts,
  double p,
  double lon,
  double lat,
  double *t) {

  double t0, t1, wt;

  /* Spatial interpolation... */
  intpol_met_space(met0, p, lon, lat, &t0);
  intpol_met_space(met1, p, lon, lat, &t1);

  /* Get weighting factor... */
  wt = (met1->time - ts) / (met1->time - met0->time);

  /* Interpolate... */
  *t = wt * (t0 - t1) + t1;
}

/*****************************************************************************/

void hamming_low_pass(
  gps_t *gps,
  double dz) {

  double ham[NZ], wsum;

  int ids, iham, iz, nham;

  /* Loop over profiles... */
  for (ids = 0; ids < gps->nds; ids++) {

    /* Calculate Hamming window coefficients... */
    nham = (int) (dz / fabs((gps->z[ids][0] - gps->z[ids][gps->nz[ids] - 1])
			    / (gps->nz[ids] - 1.0)) + 0.5);
    nham = GSL_MAX(GSL_MIN(nham, NZ), 2);
    for (iham = 0; iham < nham; iham++)
      ham[iham] = 0.54 + 0.46 * cos(M_PI * iham / (nham - 1.0));

    /* Loop over altitudes... */
    for (iz = 0; iz < gps->nz[ids]; iz++) {

      /* Initialize... */
      gps->pt[ids][iz] = ham[0] * gps->t[ids][iz];
      wsum = ham[0];

      /* Loop over filter window... */
      for (iham = 1; iham < nham; iham++) {

	/* Check array range... */
	if (iz - iham < 0 || iz + iham >= gps->nz[ids])
	  continue;

	/* Check temperature value... */
	if (!gsl_finite(gps->t[ids][iz - iham]) ||
	    !gsl_finite(gps->t[ids][iz + iham]))
	  continue;

	/* Check for tropopause... */
	if (gsl_finite(gps->th[ids]) && gps->th[ids] > 0)
	  if ((gps->z[ids][iz] >= gps->th[ids]
	       && gps->z[ids][iz - iham] < gps->th[ids])
	      || (gps->z[ids][iz] <= gps->th[ids]
		  && gps->z[ids][iz + iham] > gps->th[ids]))
	    continue;

	/* Apply Hamming filter... */
	gps->pt[ids][iz]
	  += ham[iham] * (gps->t[ids][iz - iham] + gps->t[ids][iz + iham]);
	wsum += 2 * ham[iham];
      }

      /* Calculate perturbation... */
      gps->pt[ids][iz] = gps->t[ids][iz] - gps->pt[ids][iz] / wsum;
    }
  }
}

/*****************************************************************************/

void hamming_high_pass(
  gps_t *gps,
  double dz) {

  double ham[NZ], pt[NZ], wsum;

  int ids, iham, iz, nham;

  /* Loop over profiles... */
  for (ids = 0; ids < gps->nds; ids++) {

    /* Calculate Hamming window coefficients... */
    nham = (int) (dz / fabs((gps->z[ids][0] - gps->z[ids][gps->nz[ids] - 1])
			    / (gps->nz[ids] - 1.0)) + 0.5);
    nham = GSL_MAX(GSL_MIN(nham, NZ), 2);
    for (iham = 0; iham < nham; iham++)
      ham[iham] = 0.54 + 0.46 * cos(M_PI * iham / (nham - 1.0));

    /* Loop over altitudes... */
    for (iz = 0; iz < gps->nz[ids]; iz++) {

      /* Initialize... */
      pt[iz] = ham[0] * gps->pt[ids][iz];
      wsum = ham[0];

      /* Loop over filter window... */
      for (iham = 1; iham < nham; iham++) {

	/* Check array range... */
	if (iz - iham < 0 || iz + iham >= gps->nz[ids])
	  continue;

	/* Check temperature value... */
	if (!gsl_finite(gps->t[ids][iz - iham]) ||
	    !gsl_finite(gps->t[ids][iz + iham]))
	  continue;

	/* Apply Hamming filter... */
	pt[iz]
	  += ham[iham] * (gps->pt[ids][iz - iham] + gps->pt[ids][iz + iham]);
	wsum += 2 * ham[iham];
      }

      /* Normalize... */
      pt[iz] /= wsum;
    }

    /* Set perturbation... */
    for (iz = 0; iz < gps->nz[ids]; iz++)
      gps->pt[ids][iz] = pt[iz];
  }
}

/*****************************************************************************/

void poly(
  gps_t *gps,
  int dim,
  double zmin,
  double zmax) {

  double bg[NZ];

  int ids, iz;

  /* Loop over profiles... */
  for (ids = 0; ids < gps->nds; ids++) {

    /* Set profile... */
    for (iz = 0; iz < gps->nz[ids]; iz++)
      bg[iz] = gps->pt[ids][iz];

    /* Polynomial interpolation... */
    poly_help(gps->z[ids], bg, gps->nz[ids], dim, zmin, zmax);

    /* Remove background... */
    for (iz = 0; iz < gps->nz[ids]; iz++)
      gps->pt[ids][iz] -= bg[iz];
  }
}

/*****************************************************************************/

void poly_help(
  double *xx,
  double *yy,
  int n,
  int dim,
  double xmin,
  double xmax) {

  gsl_multifit_linear_workspace *work;
  gsl_matrix *cov, *X;
  gsl_vector *c, *x, *y;

  double chisq, xx2[NZ], yy2[NZ];

  size_t i, i2, n2 = 0;

  /* Check for nan... */
  for (i = 0; i < (size_t) n; i++)
    if (xx[i] >= xmin && xx[i] <= xmax && gsl_finite(yy[i])) {
      xx2[n2] = xx[i];
      yy2[n2] = yy[i];
      n2++;
    }
  if ((int) n2 < dim) {
    for (i = 0; i < (size_t) n; i++)
      yy[i] = GSL_NAN;
    return;
  }

  /* Allocate... */
  work = gsl_multifit_linear_alloc((size_t) n2, (size_t) dim);
  cov = gsl_matrix_alloc((size_t) dim, (size_t) dim);
  X = gsl_matrix_alloc((size_t) n2, (size_t) dim);
  c = gsl_vector_alloc((size_t) dim);
  x = gsl_vector_alloc((size_t) n2);
  y = gsl_vector_alloc((size_t) n2);

  /* Compute polynomial fit... */
  for (i = 0; i < (size_t) n2; i++) {
    gsl_vector_set(x, i, xx2[i]);
    gsl_vector_set(y, i, yy2[i]);
    for (i2 = 0; i2 < (size_t) dim; i2++)
      gsl_matrix_set(X, i, i2, pow(gsl_vector_get(x, i), (double) i2));
  }
  gsl_multifit_linear(X, y, c, cov, &chisq, work);
  for (i = 0; i < (size_t) n; i++)
    yy[i] = gsl_poly_eval(c->data, (int) dim, xx[i]);

  /* Free... */
  gsl_multifit_linear_free(work);
  gsl_matrix_free(cov);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_vector_free(x);
  gsl_vector_free(y);
}

/*****************************************************************************/

void read_gps_prof(
  char *filename,
  gps_t *gps) {

  char bad[10];

  double t0, t1, zmin = 1e100, zmax = -1e100;

  int ncid, dimid, varid;

  size_t iz, nz;

  /* Open netCDF file... */
  printf("Read GPS-RO profile: %s\n", filename);
  NC(nc_open(filename, NC_NOWRITE, &ncid));

  /* Get dimensions... */
  NC(nc_inq_dimid(ncid, "MSL_alt", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &nz));
  gps->nz[gps->nds] = (int) nz;
  if (nz > NZ)
    ERRMSG("Too many altitudes!");

  /* Check data quality flag... */
  NC(nc_get_att_text(ncid, NC_GLOBAL, "bad", bad));
  if (bad[0] != '0') {
    NC(nc_close(ncid));
    return;
  }

  /* Get time... */
  NC(nc_get_att_double(ncid, NC_GLOBAL, "start_time", &t0));
  NC(nc_get_att_double(ncid, NC_GLOBAL, "stop_time", &t1));
  gps->time[gps->nds] = 0.5 * (t0 + t1) - 630720000.0;

  /* Get data... */
  NC(nc_inq_varid(ncid, "MSL_alt", &varid));
  NC(nc_get_var_double(ncid, varid, gps->z[gps->nds]));
  NC(nc_inq_varid(ncid, "Lon", &varid));
  NC(nc_get_var_double(ncid, varid, gps->lon[gps->nds]));
  NC(nc_inq_varid(ncid, "Lat", &varid));
  NC(nc_get_var_double(ncid, varid, gps->lat[gps->nds]));
  NC(nc_inq_varid(ncid, "Pres", &varid));
  NC(nc_get_var_double(ncid, varid, gps->p[gps->nds]));
  NC(nc_inq_varid(ncid, "Temp", &varid));
  NC(nc_get_var_double(ncid, varid, gps->t[gps->nds]));
  if (nc_inq_varid(ncid, "Vp", &varid) == NC_NOERR)
    NC(nc_get_var_double(ncid, varid, gps->wv[gps->nds]));

  /* Check altitude range... */
  for (iz = 0; iz < nz; iz++)
    if (gps->p[gps->nds][iz] != -999 && gps->t[gps->nds][iz] != -999) {
      zmin = GSL_MIN(zmin, gps->z[gps->nds][iz]);
      zmax = GSL_MAX(zmax, gps->z[gps->nds][iz]);
    }
  if (zmin > 5 || zmax < 35) {
    NC(nc_close(ncid));
    return;
  }

  /* Check data... */
  for (iz = 0; iz < nz; iz++)
    if (gps->lon[gps->nds][iz] == -999 ||
	gps->lat[gps->nds][iz] == -999 ||
	gps->p[gps->nds][iz] == -999 ||
	gps->t[gps->nds][iz] == -999 || gps->wv[gps->nds][iz] == -999) {
      gps->lon[gps->nds][iz] = GSL_NAN;
      gps->lat[gps->nds][iz] = GSL_NAN;
      gps->p[gps->nds][iz] = GSL_NAN;
      gps->t[gps->nds][iz] = GSL_NAN;
      gps->wv[gps->nds][iz] = GSL_NAN;
    }

  /* Convert temperature... */
  for (iz = 0; iz < nz; iz++)
    gps->t[gps->nds][iz] += 273.15;

  /* Convert water vapor... */
  for (iz = 0; iz < nz; iz++)
    gps->wv[gps->nds][iz] /= gps->p[gps->nds][iz];

  /* Close file... */
  NC(nc_close(ncid));

  /* Count profiles... */
  if ((++gps->nds) >= NDS)
    ERRMSG("Too many profiles!");
}

/*****************************************************************************/

void read_gps(
  char *filename,
  gps_t *gps) {

  int ids, ncid, dimid, varid;

  size_t start[2], count[2], nds, nz;

  /* Read netCDF file... */
  printf("Read GPS-RO file: %s\n", filename);
  NC(nc_open(filename, NC_NOWRITE, &ncid));

  /* Get dimensions... */
  NC(nc_inq_dimid(ncid, "NDS", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &nds));
  gps->nds = (int) nds;
  if (nds > NDS)
    ERRMSG("Too many profiles!");

  NC(nc_inq_dimid(ncid, "NZ", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &nz));
  if (nz > NZ)
    ERRMSG("Too many profiles!");

  /* Loop over profiles... */
  for (ids = 0; ids < gps->nds; ids++) {

    /* Set profile index... */
    start[0] = (size_t) ids;
    count[0] = 1;
    start[1] = 0;
    count[1] = nz;

    /* Set number of altitudes... */
    gps->nz[ids] = (int) nz;

    /* Read data... */
    NC(nc_inq_varid(ncid, "time", &varid));
    NC(nc_get_vara_double(ncid, varid, start, count, &gps->time[ids]));

    NC(nc_inq_varid(ncid, "z", &varid));
    NC(nc_get_vara_double(ncid, varid, start, count, gps->z[ids]));

    NC(nc_inq_varid(ncid, "lon", &varid));
    NC(nc_get_vara_double(ncid, varid, start, count, gps->lon[ids]));

    NC(nc_inq_varid(ncid, "lat", &varid));
    NC(nc_get_vara_double(ncid, varid, start, count, gps->lat[ids]));

    NC(nc_inq_varid(ncid, "p", &varid));
    NC(nc_get_vara_double(ncid, varid, start, count, gps->p[ids]));

    NC(nc_inq_varid(ncid, "t", &varid));
    NC(nc_get_vara_double(ncid, varid, start, count, gps->t[ids]));

    NC(nc_inq_varid(ncid, "wv", &varid));
    NC(nc_get_vara_double(ncid, varid, start, count, gps->wv[ids]));

    NC(nc_inq_varid(ncid, "pt", &varid));
    NC(nc_get_vara_double(ncid, varid, start, count, gps->pt[ids]));

    NC(nc_inq_varid(ncid, "th", &varid));
    NC(nc_get_vara_double(ncid, varid, start, count, &gps->th[ids]));
  }

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void read_met(
  char *filename,
  met_t *met) {

  char tstr[10];

  int ip, dimid, ncid, varid, year, mon, day, hour;

  size_t np, nx, ny;

  /* Write info... */
  printf("Read meteorological data: %s\n", filename);

  /* Get time from filename... */
  sprintf(tstr, "%.4s", &filename[strlen(filename) - 16]);
  year = atoi(tstr);
  sprintf(tstr, "%.2s", &filename[strlen(filename) - 11]);
  mon = atoi(tstr);
  sprintf(tstr, "%.2s", &filename[strlen(filename) - 8]);
  day = atoi(tstr);
  sprintf(tstr, "%.2s", &filename[strlen(filename) - 5]);
  hour = atoi(tstr);
  time2jsec(year, mon, day, hour, 0, 0, 0, &met->time);

  /* Open netCDF file... */
  NC(nc_open(filename, NC_NOWRITE, &ncid));

  /* Get dimensions... */
  NC(nc_inq_dimid(ncid, "lon", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &nx));
  if (nx > EX)
    ERRMSG("Too many longitudes!");

  NC(nc_inq_dimid(ncid, "lat", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &ny));
  if (ny > EY)
    ERRMSG("Too many latitudes!");

  NC(nc_inq_dimid(ncid, "lev", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &np));
  if (np > EP)
    ERRMSG("Too many levels!");

  /* Store dimensions... */
  met->np = (int) np;
  met->nx = (int) nx;
  met->ny = (int) ny;

  /* Get horizontal grid... */
  NC(nc_inq_varid(ncid, "lon", &varid));
  NC(nc_get_var_double(ncid, varid, met->lon));
  NC(nc_inq_varid(ncid, "lat", &varid));
  NC(nc_get_var_double(ncid, varid, met->lat));

  /* Read meteorological data... */
  read_met_help(ncid, "t", "T", met, met->t, 1.0);

  /* Read pressure levels from file... */
  NC(nc_inq_varid(ncid, "lev", &varid));
  NC(nc_get_var_double(ncid, varid, met->p));
  for (ip = 0; ip < met->np; ip++)
    met->p[ip] /= 100.;

  /* Extrapolate data for lower boundary... */
  read_met_extrapolate(met);

  /* Check ordering of pressure levels... */
  for (ip = 1; ip < met->np; ip++)
    if (met->p[ip - 1] < met->p[ip])
      ERRMSG("Pressure levels must be descending!");

  /* Create periodic boundary conditions... */
  read_met_periodic(met);

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void read_met_extrapolate(
  met_t *met) {

  int ip, ip0, ix, iy;

  /* Loop over columns... */
  for (ix = 0; ix < met->nx; ix++)
    for (iy = 0; iy < met->ny; iy++) {

      /* Find lowest valid data point... */
      for (ip0 = met->np - 1; ip0 >= 0; ip0--)
	if (!gsl_finite(met->t[ix][iy][ip0]))
	  break;

      /* Extrapolate... */
      for (ip = ip0; ip >= 0; ip--)
	met->t[ix][iy][ip] = met->t[ix][iy][ip + 1];
    }
}

/*****************************************************************************/

void read_met_help(
  int ncid,
  char *varname,
  char *varname2,
  met_t *met,
  float dest[EX][EY][EP],
  float scl) {

  static float help[EX * EY * EP];

  int ip, ix, iy, n = 0, varid;

  /* Check if variable exists... */
  if (nc_inq_varid(ncid, varname, &varid) != NC_NOERR)
    if (nc_inq_varid(ncid, varname2, &varid) != NC_NOERR)
      return;

  /* Read data... */
  NC(nc_get_var_float(ncid, varid, help));

  /* Copy and check data... */
  for (ip = 0; ip < met->np; ip++)
    for (iy = 0; iy < met->ny; iy++)
      for (ix = 0; ix < met->nx; ix++) {
	dest[ix][iy][ip] = scl * help[n++];
	if (fabs(dest[ix][iy][ip] / scl) > 1e14)
	  dest[ix][iy][ip] = GSL_NAN;
      }
}

/*****************************************************************************/

void read_met_periodic(
  met_t *met) {

  int ip, iy;

  /* Check longitudes... */
  if (!(fabs(met->lon[met->nx - 1] - met->lon[0]
	     + met->lon[1] - met->lon[0] - 360) < 0.01))
    return;

  /* Increase longitude counter... */
  if ((++met->nx) > EX)
    ERRMSG("Cannot create periodic boundary conditions!");

  /* Set longitude... */
  met->lon[met->nx - 1] = met->lon[met->nx - 2] + met->lon[1] - met->lon[0];

  /* Loop over latitudes and pressure levels... */
  for (iy = 0; iy < met->ny; iy++)
    for (ip = 0; ip < met->np; ip++)
      met->t[met->nx - 1][iy][ip] = met->t[0][iy][ip];
}

/*****************************************************************************/

void spline(
  const double *x,
  const double *y,
  const int n,
  const double *x2,
  double *y2,
  const int n2,
  const int method) {

  /* Cubic spline interpolation... */
  if (method == 1) {

    /* Allocate... */
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *s = gsl_spline_alloc(gsl_interp_cspline, (size_t) n);

    /* Interpolate profile... */
    gsl_spline_init(s, x, y, (size_t) n);
    for (int i = 0; i < n2; i++)
      if (x2[i] <= x[0])
	y2[i] = y[0];
      else if (x2[i] >= x[n - 1])
	y2[i] = y[n - 1];
      else
	y2[i] = gsl_spline_eval(s, x2[i], acc);

    /* Free... */
    gsl_spline_free(s);
    gsl_interp_accel_free(acc);
  }

  /* Linear interpolation... */
  else {
    for (int i = 0; i < n2; i++)
      if (x2[i] <= x[0])
	y2[i] = y[0];
      else if (x2[i] >= x[n - 1])
	y2[i] = y[n - 1];
      else {
	int idx = locate_irr(x, n, x2[i]);
	y2[i] = LIN(x[idx], y[idx], x[idx + 1], y[idx + 1], x2[i]);
      }
  }
}

/*****************************************************************************/

void tropopause(
  gps_t *gps) {

  double zmin;

  int ids, iz, iz2, okay;

  /* Loop over profiles... */
  for (ids = 0; ids < gps->nds; ids++) {

    /* Set default value... */
    gps->th[ids] = GSL_NAN;

    /* Set minimum altitude... */
    zmin =
      8 - 4 * fabs(cos((90 - gps->lat[ids][gps->nz[ids] / 2]) * M_PI / 180));

    /* Search tropopause (WMO definition)... */
    for (iz = 0; iz < gps->nz[ids]; iz++)
      if (gps->z[ids][iz] >= zmin && gps->z[ids][iz] <= 20.0) {
	okay = 1;
	for (iz2 = iz + 1; iz2 < gps->nz[ids]; iz2++)
	  if (gps->z[ids][iz2] - gps->z[ids][iz] <= 2.0)
	    if (!gsl_finite(gps->t[ids][iz]) ||
		!gsl_finite(gps->t[ids][iz2]) ||
		(gps->t[ids][iz2] - gps->t[ids][iz])
		/ (gps->z[ids][iz2] - gps->z[ids][iz]) < -2.0)
	      okay = 0;
	if (okay) {
	  gps->th[ids] = gps->z[ids][iz];
	  break;
	}
      }
  }
}

/*****************************************************************************/

void tropopause_spline(
  gps_t *gps,
  int met_tropo) {

  /* Loop over data sets... */
#pragma omp parallel for default(shared)
  for (int ids = 0; ids < gps->nds; ids++) {

    /* Initialize tropopause data... */
    gps->tp[ids] = NAN;
    gps->th[ids] = NAN;
    gps->tt[ids] = NAN;
    gps->tq[ids] = NAN;
    gps->tlon[ids] = NAN;
    gps->tlat[ids] = NAN;

    /* Get altitude and pressure profiles... */
    int nz = 0;
    double h[NZ], t[NZ], z[NZ], q[NZ], lon[NZ], lat[NZ];
    for (int iz = 0; iz < gps->nz[ids]; iz++)
      if (gsl_finite(gps->p[ids][iz]) && gsl_finite(gps->t[ids][iz])
	  && gsl_finite(gps->z[ids][iz])
	  && gps->z[ids][iz] >= 4.0 && gps->z[ids][iz] <= 24.0) {
	h[nz] = gps->z[ids][iz];
	t[nz] = gps->t[ids][iz];
	z[nz] = Z(gps->p[ids][iz]);
	q[nz] = gps->wv[ids][iz];
	lon[nz] = gps->lon[ids][iz];
	lat[nz] = gps->lat[ids][iz];
	if (nz > 0 && z[nz] <= z[nz - 1])
	  ERRMSG("Profiles must be ascending!");
	if ((++nz) >= NZ)
	  ERRMSG("Too many height levels!");
      }

    /* Set grid for spline interpolation... */
    double h2[200], p2[200], t2[200], z2[200], q2[200];
    for (int iz = 0; iz <= 190; iz++) {
      z2[iz] = 4.5 + 0.1 * iz;
      p2[iz] = P(z2[iz]);
    }

    /* Interpolate temperature and geopotential height profiles... */
    spline(z, t, nz, z2, t2, 191, 1);
    spline(z, h, nz, z2, h2, 191, 1);
    spline(z, q, nz, z2, q2, 191, 1);

    /* Use cold point... */
    if (met_tropo == 2) {

      /* Find minimum... */
      int iz = (int) gsl_stats_min_index(t2, 1, 171);
      if (iz > 0 && iz < 170) {
	gps->tp[ids] = p2[iz];
	gps->th[ids] = h2[iz];
	gps->tt[ids] = t2[iz];
	gps->tq[ids] = q2[iz];
      }
    }

    /* Use WMO definition... */
    else if (met_tropo == 3 || met_tropo == 4) {

      /* Find 1st tropopause... */
      int iz;
      for (iz = 0; iz <= 170; iz++) {
	int found = 1;
	for (int iz2 = iz + 1; iz2 <= iz + 20; iz2++)
	  if (LAPSE(p2[iz], t2[iz], p2[iz2], t2[iz2]) > 2.0) {
	    found = 0;
	    break;
	  }
	if (found) {
	  if (iz > 0 && iz < 170) {
	    gps->tp[ids] = p2[iz];
	    gps->th[ids] = h2[iz];
	    gps->tt[ids] = t2[iz];
	    gps->tq[ids] = q2[iz];
	  }
	  break;
	}
      }

      /* Find 2nd tropopause... */
      if (met_tropo == 4) {
	gps->tp[ids] = NAN;
	for (; iz <= 170; iz++) {
	  int found = 1;
	  for (int iz2 = iz + 1; iz2 <= iz + 10; iz2++)
	    if (LAPSE(p2[iz], t2[iz], p2[iz2], t2[iz2]) < 3.0) {
	      found = 0;
	      break;
	    }
	  if (found)
	    break;
	}
	for (; iz <= 170; iz++) {
	  int found = 1;
	  for (int iz2 = iz + 1; iz2 <= iz + 20; iz2++)
	    if (LAPSE(p2[iz], t2[iz], p2[iz2], t2[iz2]) > 2.0) {
	      found = 0;
	      break;
	    }
	  if (found) {
	    if (iz > 0 && iz < 170) {
	      gps->tp[ids] = p2[iz];
	      gps->th[ids] = h2[iz];
	      gps->tt[ids] = t2[iz];
	      gps->tq[ids] = q2[iz];
	    }
	    break;
	  }
	}
      }
    }

    else
      ERRMSG("Cannot calculate tropopause!");

    /* Find tropopause longitude and latitude... */
    if (gsl_finite(gps->th[ids]))
      for (int iz = 0; iz < nz - 1; iz++)
	if (gps->th[ids] >= h[iz] && gps->th[ids] < h[iz + 1]) {
	  gps->tlon[ids] = lon[iz];
	  gps->tlat[ids] = lat[iz];
	  break;
	}
  }
}

/*****************************************************************************/

void write_gps(
  char *filename,
  gps_t *gps) {

  static double help[NDS * NZ];

  int ids, iz, ncid, dimid[2], time_id, z_id, lon_id, lat_id, p_id, t_id,
    pt_id, wv_id, th_id, nzmax = 0;

  /* Create netCDF file... */
  printf("Write GPS-RO file: %s\n", filename);
  NC(nc_create(filename, NC_CLOBBER, &ncid));

  /* Set dimensions... */
  NC(nc_def_dim(ncid, "NDS", (size_t) gps->nds, &dimid[0]));
  for (ids = 0; ids < gps->nds; ids++)
    nzmax = GSL_MAX(nzmax, gps->nz[ids]);
  NC(nc_def_dim(ncid, "NZ", (size_t) nzmax, &dimid[1]));

  /* Add variables... */
  add_var(ncid, "time", "s", "time (seconds since 2000-01-01T00:00Z)",
	  NC_DOUBLE, dimid, &time_id, 1);
  add_var(ncid, "z", "km", "altitude", NC_FLOAT, dimid, &z_id, 2);
  add_var(ncid, "lon", "deg", "longitude", NC_FLOAT, dimid, &lon_id, 2);
  add_var(ncid, "lat", "deg", "latitude", NC_FLOAT, dimid, &lat_id, 2);
  add_var(ncid, "p", "hPa", "pressure", NC_FLOAT, dimid, &p_id, 2);
  add_var(ncid, "t", "K", "temperature", NC_FLOAT, dimid, &t_id, 2);
  add_var(ncid, "wv", "ppv", "water vapor volume mixing ratio",
	  NC_FLOAT, dimid, &wv_id, 2);
  add_var(ncid, "pt", "K", "temperature perturbation",
	  NC_FLOAT, dimid, &pt_id, 2);
  add_var(ncid, "th", "km", "tropopause height", NC_FLOAT, dimid, &th_id, 1);

  /* Leave define mode... */
  NC(nc_enddef(ncid));

  /* Write data... */
  NC(nc_put_var_double(ncid, time_id, gps->time));
  NC(nc_put_var_double(ncid, th_id, gps->th));
  for (ids = 0; ids < gps->nds; ids++)
    for (iz = 0; iz < gps->nz[ids]; iz++)
      help[ids * nzmax + iz] = gps->z[ids][iz];
  NC(nc_put_var_double(ncid, z_id, help));
  for (ids = 0; ids < gps->nds; ids++)
    for (iz = 0; iz < gps->nz[ids]; iz++)
      help[ids * nzmax + iz] = gps->lon[ids][iz];
  NC(nc_put_var_double(ncid, lon_id, help));
  for (ids = 0; ids < gps->nds; ids++)
    for (iz = 0; iz < gps->nz[ids]; iz++)
      help[ids * nzmax + iz] = gps->lat[ids][iz];
  NC(nc_put_var_double(ncid, lat_id, help));
  for (ids = 0; ids < gps->nds; ids++)
    for (iz = 0; iz < gps->nz[ids]; iz++)
      help[ids * nzmax + iz] = gps->p[ids][iz];
  NC(nc_put_var_double(ncid, p_id, help));
  for (ids = 0; ids < gps->nds; ids++)
    for (iz = 0; iz < gps->nz[ids]; iz++)
      help[ids * nzmax + iz] = gps->t[ids][iz];
  NC(nc_put_var_double(ncid, t_id, help));
  for (ids = 0; ids < gps->nds; ids++)
    for (iz = 0; iz < gps->nz[ids]; iz++)
      help[ids * nzmax + iz] = gps->wv[ids][iz];
  NC(nc_put_var_double(ncid, wv_id, help));
  for (ids = 0; ids < gps->nds; ids++)
    for (iz = 0; iz < gps->nz[ids]; iz++)
      help[ids * nzmax + iz] = gps->pt[ids][iz];
  NC(nc_put_var_double(ncid, pt_id, help));

  /* Close file... */
  NC(nc_close(ncid));
}
