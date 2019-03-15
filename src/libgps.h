#include <netcdf.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_poly.h>
#include "jurassic.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of pressure levels for meteorological data. */
#define EP 73

/*! Maximum number of longitudes for meteorological data. */
#define EX 721

/*! Maximum number of latitudes for meteorological data. */
#define EY 361

/*! Maximum number of GPS-RO profiles. */
#define NDS 10000

/*! Maximum number of altitudes per GPS-RO profile. */
#define NZ 5000

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/*! Execute netCDF library command and check result. */
#define NC(cmd) {				     \
    if((cmd)!=NC_NOERR)				     \
      ERRMSG(nc_strerror(cmd));			     \
  }

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/*! GPS-RO profile data. */
typedef struct {

  /*! Number of profiles. */
  int nds;

  /*! Number of altitudes per profile. */
  int nz[NDS];

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time[NDS];

  /*! Altitude [km]. */
  double z[NDS][NZ];

  /*! Longitude [deg]. */
  double lon[NDS][NZ];

  /*! Latitude [deg]. */
  double lat[NDS][NZ];

  /*! Pressure [hPa]. */
  double p[NDS][NZ];

  /*! Temperature [K]. */
  double t[NDS][NZ];

  /*! Water vapor volume mixing ratio [ppm]. */
  double wv[NDS][NZ];

  /*! Temperature perturbation [K]. */
  double pt[NDS][NZ];

  /*! Tropopause height [km]. */
  double th[NDS];

} gps_t;

/*! Meteorological data. */
typedef struct {

  /*! Time [s]. */
  double time;

  /*! Number of longitudes. */
  int nx;

  /*! Number of latitudes. */
  int ny;

  /*! Number of pressure levels. */
  int np;

  /*! Longitude [deg]. */
  double lon[EX];

  /*! Latitude [deg]. */
  double lat[EY];

  /*! Pressure [hPa]. */
  double p[EP];

  /*! Temperature [K]. */
  float t[EX][EY][EP];

} met_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Add variable to netCDF file. */
void add_var(
  int ncid,
  const char *varname,
  const char *unit,
  const char *longname,
  int type,
  int dimid[],
  int *varid,
  int ndims);

/*! Detrending by means of meteo data. */
void detrend_met(
  gps_t * gps,
  char *metbase,
  double dt_met);

/*! Calculate horizontal Gaussian mean to extract perturbations. */
void gauss(
  gps_t * gps,
  double dx,
  double dy);

/*! Interpolate GPS data to regular altitude grid. */
void grid_gps(
  gps_t * gps,
  double zmin,
  double zmax,
  int nz);

/*! Get meteorological data for given timestep. */
void get_met(
  char *metbase,
  double dt_met,
  double t,
  met_t * met0,
  met_t * met1);

/*! Get meteorological data for timestep. */
void get_met_help(
  double t,
  int direct,
  char *metbase,
  double dt_met,
  char *filename);

/*! Linear interpolation of 3-D meteorological data. */
void intpol_met_3d(
  float array[EX][EY][EP],
  int ip,
  int ix,
  int iy,
  double wp,
  double wx,
  double wy,
  double *var);

/*! Spatial interpolation of meteorological data. */
void intpol_met_space(
  met_t * met,
  double p,
  double lon,
  double lat,
  double *t);

/*! Temporal interpolation of meteorological data. */
void intpol_met_time(
  met_t * met0,
  met_t * met1,
  double ts,
  double p,
  double lon,
  double lat,
  double *t);

/*! Apply vertical Hamming filter to extract perturbations. */
void hamming_low_pass(
  gps_t * gps,
  double dz);

/*! Apply vertical Hamming filter to reduce noise. */
void hamming_high_pass(
  gps_t * gps,
  double dz);

/*! Remove polynomial fit from perturbation profile. */
void poly(
  gps_t * gps,
  int dim,
  double zmin,
  double zmax);

/*! Auxiliary function for polynomial interpolation. */
void poly_help(
  double *xx,
  double *yy,
  int n,
  int dim,
  double xmin,
  double xmax);

/*! Read GPS-RO profile. */
void read_gps_prof(
  char *filename,
  gps_t * gps);

/*! Read GPS-RO data file. */
void read_gps(
  char *filename,
  gps_t * gps);

/*! Read meteorological data file. */
void read_met(
  char *filename,
  met_t * met);

/*! Extrapolate meteorological data at lower boundary. */
void read_met_extrapolate(
  met_t * met);

/*! Read and convert variable from meteorological data file. */
void read_met_help(
  int ncid,
  char *varname,
  char *varname2,
  met_t * met,
  float dest[EX][EY][EP],
  float scl);

/*! Create meteorological data with periodic boundary conditions. */
void read_met_periodic(
  met_t * met);

/*! Find tropopause height. */
void tropopause(
  gps_t * gps);

/*! Write GPS-RO data file. */
void write_gps(
  char *filename,
  gps_t * gps);
