/******************************************************************************
 * fixproj is designed to take spherical coordinates from a NetCDF file       *
 * (lon,lat), transform the coordinates with a PROJ.4 cartographic projection *
 * reference, and save the results (x,y). By default fixproj will attempt to  *
 * read the projection reference from the NetCDF file's global                *
 * CoordinateProjection attribute, but if a reference is not found a new one  *
 * will be generated to optimally fit the bounds of the data.                 *
 *                                                                            *
 * Available command line options are:                                        *
 *  -p proj                                                                   *
 *      Provide a geographic projection reference string to be used.          *
 *                                                                            *
 *  -i                                                                        *
 *      Perform an inverse projection.                                        *
 *                                                                            *
 *  -v longitude_var latitude_var x_var y_var                                 *
 *      Provide alternative names for the coordinate variables.               *
 *                                                                            *
 *  -g                                                                        *
 *      Force the generation of a new projection reference.                   *
 *                                                                            *
 *  -s                                                                        *
 *      Suppresses the write-back of the used projection reference.           *
 *                                                                            *
 * Author: Tristan.Losier@unb.ca                                              *
 ******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <netcdf.h>
#include <proj_api.h>

/******************************************************************************
 *  Header Data                                                               *
 ******************************************************************************/
#define USAGE_STRING  "Command usage:\n    fixproj [-p projection] [-i] [-v lon lat x y] [-g] [-s] file_to_process.nc"
#define PROJ_ATT_NAME "CoordinateProjection"
#define PROJ_FORMAT   "+proj=lcc +lon_0=%0.2f +lat_0=%0.2f +lat_1=%0.1f +lat_2=%0.1f +ellps=WGS84"

#define ARG_PROJ 'p'
#define ARG_INV  'i'
#define ARG_VAR  'v'
#define ARG_GEN  'g'
#define ARG_SUPR 's'

#define ERR_MEM    0
#define ERR_NETCDF 1
#define ERR_PROJ   2

projPJ proj_src;
projPJ proj_dst;
int    proj_inv;

typedef struct _ncdata {
    int    ncid;
    int    lonid;
    int    latid;
    int    xid;
    int    yid;
    size_t count;
} ncdata;

void save_proj_ref(int ncid, const char * proj);
ncdata * init_data(const char * file,
                   const char * lon,
                   const char * lat,
                   const char * x,
                   const char * y);
char * get_proj_ref(ncdata * data, int gen);
void init_proj(char * proj);
void project(ncdata * data);
void sanitize_angles(int count, ...);
void handle_error(int source, int status);

/******************************************************************************
 *  Program entry point                                                       *
 ******************************************************************************/
int main(int argc, char ** argv)
{
    char * lonvar = "lon";
    char * latvar = "lat";
    char * xvar = "x";
    char * yvar = "y";
    char * ncdf = NULL;
    char * proj = NULL;
    int gen_proj = 0;
    int writeback = 1;
    int i;
    ncdata * data;

    // Read input parameters
    for (i = 0; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
            case ARG_PROJ:
                if (argc - i < 2) break;
                proj = argv[++i];
                break;

            case ARG_INV:
                proj_inv = 1;
                break;

            case ARG_VAR:
                if (argc - i < 5) break;
                lonvar = argv[++i];
                latvar = argv[++i];
                xvar = argv[++i];
                yvar = argv[++i];
                break;

            case ARG_GEN:
                gen_proj = 1;
                break;

            case ARG_SUPR:
                writeback = 0;
            }
        }
        else if (i > 0) ncdf = argv[i];
    }

    // Make sure the user provided necessary arguments
    if (ncdf == NULL) {
        fprintf(stderr, "ERROR: NetCDF input file not specified\n");
        fprintf(stderr, "%s\n", USAGE_STRING);
        exit(EXIT_FAILURE);
    }

    // Do the projection
    data = init_data(ncdf, lonvar, latvar, xvar, yvar);
    if (proj == NULL) proj = get_proj_ref(data, gen_proj);
    if (writeback) save_proj_ref(data->ncid, proj);
    init_proj(proj);
    project(data);

    // Done, close the NetCDF file and exit
    pj_free(proj_src);
    pj_free(proj_dst);
    nc_close(data->ncid);
    exit(EXIT_SUCCESS);
}

/******************************************************************************
 *  Save a projection reference string to the provided NetCDF file            *
 ******************************************************************************/
void save_proj_ref(int ncid, const char * proj)
{
    int status = nc_redef(ncid);
    handle_error(ERR_NETCDF, status);
    status = nc_put_att_text(ncid, NC_GLOBAL, PROJ_ATT_NAME, strlen(proj), proj);
    handle_error(ERR_NETCDF, status);
    status = nc_enddef(ncid);
    handle_error(ERR_NETCDF, status);
}

/******************************************************************************
 *  Create and initialize a ncdata object using the NetCDF file provided      *
 ******************************************************************************/
ncdata * init_data(const char * file,
                   const char * lon,
                   const char * lat,
                   const char * x,
                   const char * y)
{
    int status;
    ncdata * data = malloc(sizeof(ncdata));
    if (data == NULL) handle_error(ERR_MEM, 0);

    // Open NetCDF file
    status = nc_open(file, NC_WRITE, &data->ncid);
    handle_error(ERR_NETCDF, status);

    // Get NetCDF variable ids for the coordinate variables
    status = nc_inq_varid(data->ncid, lon, &data->lonid);
    handle_error(ERR_NETCDF, status);
    status = nc_inq_varid(data->ncid, lat, &data->latid);
    handle_error(ERR_NETCDF, status);
    status = nc_inq_varid(data->ncid, x, &data->xid);
    handle_error(ERR_NETCDF, status);
    status = nc_inq_varid(data->ncid, y, &data->yid);
    handle_error(ERR_NETCDF, status);

    // Find the number of points being transformed
    int ndim;
    status = nc_inq_varndims(data->ncid, data->xid, &ndim);
    handle_error(ERR_NETCDF, status);
    int dimids[ndim];
    status = nc_inq_vardimid(data->ncid, data->xid, dimids);
    handle_error(ERR_NETCDF, status);
    data->count = 1;
    int i;
    for (i = 0; i < ndim; i++) {
        size_t dimlen;
        status = nc_inq_dimlen(data->ncid, dimids[i], &dimlen);
        handle_error(ERR_NETCDF, status);
        data->count *= dimlen;
    }

    return data;
}

/******************************************************************************
 *  Get a projection reference from the provided NetCDF file. Either from     *
 *  the PROJ_ATT_NAME global attribute, or by generating one to fit the       *
 *  boundaries of the data contained within the file.                         *
 ******************************************************************************/
char * get_proj_ref(ncdata * data, int gen)
{
    char * proj;

    // If the user did not provide a projection reference via the command line
    // we first check the NetCDF file's CoordinateProjection global attribute
    size_t length;
    int status = nc_inq_attlen(data->ncid, NC_GLOBAL, PROJ_ATT_NAME, &length);
    if (status == NC_NOERR && length > 0 && !gen) {
        printf("Using the %s global attribute.\n", PROJ_ATT_NAME);
        proj = malloc(length + 1); // +1 for trailing NULL
        if (proj == NULL) handle_error(ERR_MEM, 0);
        status = nc_get_att_text(data->ncid, NC_GLOBAL, PROJ_ATT_NAME, proj);
        handle_error(ERR_NETCDF, status);
        proj[length] = '\0';
        return proj;
    }

    // If the NetCDF file does not have a projection reference we create one
    // that produces minimal distortion across the mesh domain
    printf("Generating new projection reference.\n");
    float * lon = malloc(data->count*sizeof(float));
    float * lat = malloc(data->count*sizeof(float));
    if (lon == NULL || lat == NULL) handle_error(ERR_MEM, 0);
    status = nc_get_var_float(data->ncid, data->lonid, lon);
    handle_error(ERR_NETCDF, status);
    status = nc_get_var_float(data->ncid, data->latid, lat);
    handle_error(ERR_NETCDF, status);

    double minlon = lon[0];
    double maxlon = lon[0];
    double minlat = lat[0];
    double maxlat = lat[0];
    int i;
    for (i = 0; i < data->count; i++) {
        if (lon[i] < minlon) minlon = lon[i];
        if (lon[i] > maxlon) maxlon = lon[i];
        if (lat[i] < minlat) minlat = lat[i];
        if (lat[i] > maxlat) maxlat = lat[i];
    }
    sanitize_angles(4, &minlon, &maxlon, &minlat, &maxlat);
    double xc = (maxlon - minlon)/2.0 + minlon;
    double yc = (maxlat - minlat)/2.0 + minlat;
    double y1 = yc - (maxlat - minlat)*2.0/6.0;
    double y2 = yc + (maxlat - minlat)*2.0/6.0;

    char tempstr[90];
    sprintf(tempstr, PROJ_FORMAT, xc, yc, y1, y2);
    proj = malloc(strlen(tempstr) + 1);
    if (proj == NULL) handle_error(ERR_MEM, 0);
    strcpy(proj, tempstr);

    free(lon);
    free(lat);
    return proj;
}

/******************************************************************************
 *  Initialize the global forward and reverse projection references           *
 ******************************************************************************/
void init_proj(char * proj)
{
    // Initialization routine for the projection lib
    proj_dst = pj_init_plus(proj);
    if (proj_dst == NULL) handle_error(ERR_PROJ, *pj_get_errno_ref());
    proj_src = pj_latlong_from_proj(proj_dst);
    if (proj_src == NULL) handle_error(ERR_PROJ, *pj_get_errno_ref());

    // Swap projection references if we need a reverse projection
    if (proj_inv) {
        projPJ tmp = proj_src;
        proj_src = proj_dst;
        proj_dst = tmp;
    }

    printf("Source projection:      %s\n", pj_get_def(proj_src, 0));
    printf("Destination projection: %s\n", pj_get_def(proj_dst, 0));
}

/******************************************************************************
 *  Use the global projection references to transform the data in the         *
 *  NetCDF file represented by the provided ncdata reference                  *
 ******************************************************************************/
void project(ncdata * data)
{
    int status, i;

    // Allocate memory for the transformation
    float  * nc_x = malloc(data->count*sizeof(float));
    float  * nc_y = malloc(data->count*sizeof(float));
    double * pj_x = malloc(data->count*sizeof(double));
    double * pj_y = malloc(data->count*sizeof(double));
    if (nc_x == NULL || nc_y == NULL || pj_x == NULL || pj_y == NULL)
        handle_error(ERR_MEM, 0);

    // Get the source coordinates
    if (proj_inv) {
        status = nc_get_var_float(data->ncid, data->xid, nc_x);
        handle_error(ERR_NETCDF, status);
        status = nc_get_var_float(data->ncid, data->yid, nc_y);
        handle_error(ERR_NETCDF, status);
    } else {
        status = nc_get_var_float(data->ncid, data->lonid, nc_x);
        handle_error(ERR_NETCDF, status);
        status = nc_get_var_float(data->ncid, data->latid, nc_y);
        handle_error(ERR_NETCDF, status);
    }

    for (i = 0; i < data->count; i++) {
        pj_x[i] = (double)nc_x[i];
        pj_y[i] = (double)nc_y[i];
        if (!proj_inv) {
            pj_x[i] *= DEG_TO_RAD;
            pj_y[i] *= DEG_TO_RAD;
        }
    }

    // Do the transformation
    status = pj_transform(proj_src, proj_dst, data->count, 1, pj_x, pj_y, NULL);
    handle_error(ERR_PROJ, status);
    printf("Done transforming %d points\n", (int)data->count);

    // Save the results
    for (i = 0; i < data->count; i++) {
        if (proj_inv) {
            pj_x[i] *= RAD_TO_DEG;
            pj_y[i] *= RAD_TO_DEG;
        }
        nc_x[i] = (float)pj_x[i];
        nc_y[i] = (float)pj_y[i];
    }

    if (proj_inv) {
        status = nc_put_var_float(data->ncid, data->lonid, nc_x);
        handle_error(ERR_NETCDF, status);
        status = nc_put_var_float(data->ncid, data->latid, nc_y);
        handle_error(ERR_NETCDF, status);
    } else {
        status = nc_put_var_float(data->ncid, data->xid, nc_x);
        handle_error(ERR_NETCDF, status);
        status = nc_put_var_float(data->ncid, data->yid, nc_y);
        handle_error(ERR_NETCDF, status);
    }

    free(nc_x);
    free(nc_y);
    free(pj_x);
    free(pj_y);
}

/******************************************************************************
 *  Take a list of angles in degrees (double *), and adjusts them to fall     *
 *  within the range (-180, 180)                                              *
 ******************************************************************************/
void sanitize_angles(int count, ...)
{
    int i;
    va_list vars;
    va_start(vars, count);
    for (i = 0; i < count; i++) {
        double * val = va_arg(vars, double *);
        while (*val >  180.0) *val -= 360.0;
        while (*val < -180.0) *val += 360.0;
    }
    va_end(vars);
}

/******************************************************************************
 *  Check for standard errors, and print descriptive error text if needed     *
 ******************************************************************************/
void handle_error(int source, int status)
{
    switch(source) {

    case ERR_MEM:
        fprintf(stderr, "ERROR: Failed to allocate memory\n");
        exit(EXIT_FAILURE);
        break;

    case ERR_NETCDF:
        if (status != NC_NOERR) {
            fprintf(stderr, "%s\n", nc_strerror(status));
            exit(EXIT_FAILURE);
        }
        break;

    case ERR_PROJ:
        if (status != 0) {
            fprintf(stderr, "%s\n", pj_strerrno(status));
            exit(EXIT_FAILURE);
        }
        break;
    }
}

