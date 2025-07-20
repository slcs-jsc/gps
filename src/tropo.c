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
  Extract tropopause data.
*/

#include "libgps.h"

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  gps_t *gps;
  
  double clp_tlon[NDS], clp_tlat[NDS], wmo_1st_tlon[NDS], wmo_1st_tlat[NDS],
    wmo_2nd_tlon[NDS], wmo_2nd_tlat[NDS];

  float clp_th[NDS], clp_tp[NDS], clp_tq[NDS], clp_tt[NDS], wmo_1st_th[NDS],
    wmo_1st_tp[NDS], wmo_1st_tq[NDS], wmo_1st_tt[NDS], wmo_2nd_th[NDS],
    wmo_2nd_tp[NDS], wmo_2nd_tq[NDS], wmo_2nd_tt[NDS];

  /* Allocate... */
  ALLOC(gps, gps_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <out.nc> <gps1.nc> [<gps2.nc> ...]");

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

  /* ------------------------------------------------------------
     Get tropopause data...
     ------------------------------------------------------------ */

  /* Cold point... */
  tropopause_spline(gps, 2);
  for (int ids = 0; ids < gps->nds; ids++) {
    clp_tlon[ids] = gps->tlon[ids];
    clp_tlat[ids] = gps->tlat[ids];
    clp_tp[ids] = (float) gps->tp[ids];
    clp_th[ids] = (float) gps->th[ids];
    clp_tt[ids] = (float) gps->tt[ids];
    clp_tq[ids] = (float) gps->tq[ids];
  }

  /* WMO 1st tropopause... */
  tropopause_spline(gps, 3);
  for (int ids = 0; ids < gps->nds; ids++) {
    wmo_1st_tlon[ids] = gps->tlon[ids];
    wmo_1st_tlat[ids] = gps->tlat[ids];
    wmo_1st_tp[ids] = (float) gps->tp[ids];
    wmo_1st_th[ids] = (float) gps->th[ids];
    wmo_1st_tt[ids] = (float) gps->tt[ids];
    wmo_1st_tq[ids] = (float) gps->tq[ids];
  }

  /* WMO 2nd tropopause... */
  tropopause_spline(gps, 4);
  for (int ids = 0; ids < gps->nds; ids++) {
    wmo_2nd_tlon[ids] = gps->tlon[ids];
    wmo_2nd_tlat[ids] = gps->tlat[ids];
    wmo_2nd_tp[ids] = (float) gps->tp[ids];
    wmo_2nd_th[ids] = (float) gps->th[ids];
    wmo_2nd_tt[ids] = (float) gps->tt[ids];
    wmo_2nd_tq[ids] = (float) gps->tq[ids];
  }

  /* ------------------------------------------------------------
     Write netCDF files...
     ------------------------------------------------------------ */

  /* Create netCDF file... */
  int ncid, varid, dimid[10];
  printf("Write tropopause file: %s\n", argv[2]);
  NC(nc_create(argv[2], NC_CLOBBER, &ncid));

  /* Set dimensions... */
  NC(nc_def_dim(ncid, "NDS", (size_t) gps->nds, &dimid[0]));

  /* Add variables... */
  add_var(ncid, "time", "s", "time (seconds since 2000-01-01T00:00Z)",
	  NC_DOUBLE, dimid, &varid, 1);

  add_var(ncid, "clp_lat", "degrees_north", "cold point latitude", NC_DOUBLE,
	  dimid, &varid, 1);
  add_var(ncid, "clp_lon", "degrees_east", "cold point longitude", NC_DOUBLE,
	  dimid, &varid, 1);
  add_var(ncid, "clp_z", "km", "cold point height", NC_FLOAT, dimid, &varid,
	  1);
  add_var(ncid, "clp_p", "hPa", "cold point pressure", NC_FLOAT, dimid,
	  &varid, 1);
  add_var(ncid, "clp_t", "K", "cold point temperature", NC_FLOAT, dimid,
	  &varid, 1);
  add_var(ncid, "clp_q", "ppv", "cold point water vapor", NC_FLOAT, dimid,
	  &varid, 1);

  add_var(ncid, "wmo_1st_lat", "degrees_north", "WMO 1st tropopause latitude",
	  NC_DOUBLE, dimid, &varid, 1);
  add_var(ncid, "wmo_1st_lon", "degrees_east", "WMO 1st tropopause longitude",
	  NC_DOUBLE, dimid, &varid, 1);
  add_var(ncid, "wmo_1st_z", "km", "WMO 1st tropopause height", NC_FLOAT,
	  dimid, &varid, 1);
  add_var(ncid, "wmo_1st_p", "hPa", "WMO 1st tropopause pressure", NC_FLOAT,
	  dimid, &varid, 1);
  add_var(ncid, "wmo_1st_t", "K", "WMO 1st tropopause temperature", NC_FLOAT,
	  dimid, &varid, 1);
  add_var(ncid, "wmo_1st_q", "ppv", "WMO 1st tropopause water vapor",
	  NC_FLOAT, dimid, &varid, 1);

  add_var(ncid, "wmo_2nd_lat", "degrees_north", "WMO 2nd tropopause latitude",
	  NC_DOUBLE, dimid, &varid, 1);
  add_var(ncid, "wmo_2nd_lon", "degrees_east", "WMO 2nd tropopause longitude",
	  NC_DOUBLE, dimid, &varid, 1);
  add_var(ncid, "wmo_2nd_z", "km", "WMO 2nd tropopause height", NC_FLOAT,
	  dimid, &varid, 1);
  add_var(ncid, "wmo_2nd_p", "hPa", "WMO 2nd tropopause pressure", NC_FLOAT,
	  dimid, &varid, 1);
  add_var(ncid, "wmo_2nd_t", "K", "WMO 2nd tropopause temperature", NC_FLOAT,
	  dimid, &varid, 1);
  add_var(ncid, "wmo_2nd_q", "ppv", "WMO 2nd tropopause water vapor",
	  NC_FLOAT, dimid, &varid, 1);

  /* Leave define mode... */
  NC(nc_enddef(ncid));

  /* Write data... */
  NC(nc_inq_varid(ncid, "time", &varid));
  NC(nc_put_var_double(ncid, varid, gps->time));

  NC(nc_inq_varid(ncid, "clp_lat", &varid));
  NC(nc_put_var_double(ncid, varid, clp_tlat));
  NC(nc_inq_varid(ncid, "clp_lon", &varid));
  NC(nc_put_var_double(ncid, varid, clp_tlon));
  NC(nc_inq_varid(ncid, "clp_z", &varid));
  NC(nc_put_var_float(ncid, varid, clp_th));
  NC(nc_inq_varid(ncid, "clp_p", &varid));
  NC(nc_put_var_float(ncid, varid, clp_tp));
  NC(nc_inq_varid(ncid, "clp_t", &varid));
  NC(nc_put_var_float(ncid, varid, clp_tt));
  NC(nc_inq_varid(ncid, "clp_q", &varid));
  NC(nc_put_var_float(ncid, varid, clp_tq));

  NC(nc_inq_varid(ncid, "wmo_1st_lat", &varid));
  NC(nc_put_var_double(ncid, varid, wmo_1st_tlat));
  NC(nc_inq_varid(ncid, "wmo_1st_lon", &varid));
  NC(nc_put_var_double(ncid, varid, wmo_1st_tlon));
  NC(nc_inq_varid(ncid, "wmo_1st_z", &varid));
  NC(nc_put_var_float(ncid, varid, wmo_1st_th));
  NC(nc_inq_varid(ncid, "wmo_1st_p", &varid));
  NC(nc_put_var_float(ncid, varid, wmo_1st_tp));
  NC(nc_inq_varid(ncid, "wmo_1st_t", &varid));
  NC(nc_put_var_float(ncid, varid, wmo_1st_tt));
  NC(nc_inq_varid(ncid, "wmo_1st_q", &varid));
  NC(nc_put_var_float(ncid, varid, wmo_1st_tq));

  NC(nc_inq_varid(ncid, "wmo_2nd_lat", &varid));
  NC(nc_put_var_double(ncid, varid, wmo_2nd_tlat));
  NC(nc_inq_varid(ncid, "wmo_2nd_lon", &varid));
  NC(nc_put_var_double(ncid, varid, wmo_2nd_tlon));
  NC(nc_inq_varid(ncid, "wmo_2nd_z", &varid));
  NC(nc_put_var_float(ncid, varid, wmo_2nd_th));
  NC(nc_inq_varid(ncid, "wmo_2nd_p", &varid));
  NC(nc_put_var_float(ncid, varid, wmo_2nd_tp));
  NC(nc_inq_varid(ncid, "wmo_2nd_t", &varid));
  NC(nc_put_var_float(ncid, varid, wmo_2nd_tt));
  NC(nc_inq_varid(ncid, "wmo_2nd_q", &varid));
  NC(nc_put_var_float(ncid, varid, wmo_2nd_tq));

  /* Close file... */
  NC(nc_close(ncid));

  /* Free... */
  free(gps);

  return EXIT_SUCCESS;
}
