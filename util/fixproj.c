/******************************************************************************
 * This program applies a cartographic projection to the spherical coordinate *
 * values contained in the provided NetCDF file. If the projection reference  *
 * was not provided via the command line with the -p option, this program     *
 * will attempt to read it from the NetCDF file's "CoordinateProjection"      *
 * global attribute.                                                          *
 *                                                                            *
 * The -i argument causes the program to perform an inverse projection.       *
 *                                                                            *
 * The -v argument provides alternative names for the Cartesian and spherical *
 * coordinate variables, in the order:                                        *
 *   -v longitude_var latitude_var x_var y_var                                *
 *                                                                            *
 * Author: Tristan.Losier@unb.ca                                              *
 ******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <proj_api.h>

#define USAGE_STRING "Command usage:\n    fixproj [-p projection] [-i] [-v lon lat x y] file_to_process.nc"
#define PROJ_ATT_NAME "CoordinateProjection"

#define ARG_PROJ 'p'
#define ARG_INV  'i'
#define ARG_VAR  'v'

#define ERR_MEM    0
#define ERR_NETCDF 1
#define ERR_PROJ   2

projPJ proj_src;
projPJ proj_dst;
int proj_inv;

void init_proj(char * proj, int ncid);
void project(int ncid, int inx, int iny, int outx, int outy);
void handle_error(int source, int status);

int main(int argc, char ** argv)
{
    char * lonvar = "lon";
    char * latvar = "lat";
    char * xvar = "x";
    char * yvar = "y";
    char * ncdf = NULL;
    char * proj = NULL;
    int ncid, status, i;
    int lonid, latid, xid, yid;

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
            }
        }
        else if (i > 0) ncdf = argv[i];
    }

    // Open NetCDF file
    if (ncdf == NULL) {
        fprintf(stderr, "ERROR: NetCDF input file not specified\n");
        printf("%s\n", USAGE_STRING);
        exit(EXIT_FAILURE);
    }
    status = nc_open(ncdf, NC_WRITE, &ncid);
    handle_error(ERR_NETCDF, status);

    // Get NetCDF variable ids for the coordinate variables
    status = nc_inq_varid(ncid, lonvar, &lonid);
    handle_error(ERR_NETCDF, status);
    status = nc_inq_varid(ncid, latvar, &latid);
    handle_error(ERR_NETCDF, status);
    status = nc_inq_varid(ncid, xvar, &xid);
    handle_error(ERR_NETCDF, status);
    status = nc_inq_varid(ncid, yvar, &yid);
    handle_error(ERR_NETCDF, status);

    // Do the projection
    init_proj(proj, ncid);
    if (!proj_inv)
        project(ncid, lonid, latid, xid, yid);
    else
        project(ncid, xid, yid, lonid, latid);

    // Done, close the NetCDF file and exit
    pj_free(proj_src);
    pj_free(proj_dst);
    status = nc_close(ncid);
    handle_error(ERR_NETCDF, status);
    exit(EXIT_SUCCESS);
}

void init_proj(char * proj, int ncid)
{
    int status;
    if (proj == NULL) {
        // If the user did not provide a projection reference via the command line
        // we read it from the NetCDF file's CoordinateProjection global attribute
        printf("No projection reference specified, using the %s global attribute\n", PROJ_ATT_NAME);
        size_t length;
        status = nc_inq_attlen(ncid, NC_GLOBAL, PROJ_ATT_NAME, &length);
        handle_error(ERR_NETCDF, status);

        proj = malloc(length + 1); // +1 for trailing NULL
        if (proj == NULL) handle_error(ERR_MEM, 0);
        status = nc_get_att_text(ncid, NC_GLOBAL, PROJ_ATT_NAME, proj);
        handle_error(ERR_NETCDF, status);

        proj[length] = '\0';
    } else {
        // If they did, update the file's projection reference attribute instead
        status = nc_redef(ncid);
        handle_error(ERR_NETCDF, status);
        status = nc_put_att_text(ncid, NC_GLOBAL, PROJ_ATT_NAME, strlen(proj), proj);
        handle_error(ERR_NETCDF, status);
        status = nc_enddef(ncid);
        handle_error(ERR_NETCDF, status);
    }

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

void project(int ncid, int inx, int iny, int outx, int outy)
{
    int status, i;
    size_t count = 1;

    // Find the number of points being transformed
    int ndim;
    status = nc_inq_varndims(ncid, inx, &ndim);
    handle_error(ERR_NETCDF, status);
    int dimids[ndim];
    status = nc_inq_vardimid(ncid, inx, dimids);
    handle_error(ERR_NETCDF, status);
    for (i = 0; i < ndim; i++) {
        size_t dimlen;
        status = nc_inq_dimlen(ncid, dimids[i], &dimlen);
        handle_error(ERR_NETCDF, status);
        count *= dimlen;
    }

    // Allocate memory for the transformation
    float * nc_x = malloc(count*sizeof(float));
    float * nc_y = malloc(count*sizeof(float));
    double * pj_x = malloc(count*sizeof(double));
    double * pj_y = malloc(count*sizeof(double));
    if (nc_x == NULL || nc_y == NULL || pj_x == NULL || pj_y == NULL)
        handle_error(ERR_MEM, 0);

    // Get the source coordinates
    status = nc_get_var_float(ncid, inx, nc_x);
    handle_error(ERR_NETCDF, status);
    status = nc_get_var_float(ncid, iny, nc_y);
    handle_error(ERR_NETCDF, status);
    for (i = 0; i < count; i++) {
        pj_x[i] = (double)nc_x[i];
        pj_y[i] = (double)nc_y[i];
        if (!proj_inv) {
            pj_x[i] *= DEG_TO_RAD;
            pj_y[i] *= DEG_TO_RAD;
        }
    }

    // Do the transformation
    status = pj_transform(proj_src, proj_dst, count, 1, pj_x, pj_y, NULL);
    handle_error(ERR_PROJ, status);
    printf("Done transforming %d points\n", (int)count);

    // Save the results
    for (i = 0; i < count; i++) {
        if (proj_inv) {
            pj_x[i] *= RAD_TO_DEG;
            pj_y[i] *= RAD_TO_DEG;
        }
        nc_x[i] = (float)pj_x[i];
        nc_y[i] = (float)pj_y[i];
    }
    status = nc_put_var_float(ncid, outx, nc_x);
    handle_error(ERR_NETCDF, status);
    status = nc_put_var_float(ncid, outy, nc_y);
    handle_error(ERR_NETCDF, status);

    free(nc_x);
    free(nc_y);
    free(pj_x);
    free(pj_y);
}

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

