/******************************************************************************
 * This program converts a particle_release.dat file to NetCDF format,        *
 * suitable for use by PTrack.                                                *
 *                                                                            *
 * The -p argument is required, and is used to provide a PROJ.4 compatible    *
 * geographic projection reference.                                           *
 *                                                                            *
 * The -o option may be used to specify an output file name. By default the   *
 * file will have the same name as the input file, except with a ".nc"        *
 * extension instead of ".dat".                                               *
 *                                                                            *
 * The -t option tells the program to read start and start dates formatted    *
 * as "dd/mm/yyyy hh:mm:ss" instead of a Modified Julian Date.                *
 ******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <proj_api.h>
#include <netcdf.h>

#define USAGE_STRING "Command usage:\n    genseed [-o outfile.nc] [-t] -p projection particle_release.dat"

#define ARG_PROJ 'p'
#define ARG_OUT  'o'
#define ARG_TIME 't'

typedef struct _point {
    struct _point * next;
    float x;
    float y;
    float z;
    float start;
    float end;
    int id;
} point;

void project(const char * proj, point * plist);
void create_outfile(const char * file,
                           int * pfid,
                           int * pnumid,
                           int * pxid,
                           int * pyid,
                           int * pzid,
                           int * prelid,
                           int * pendid);
point * read_file(const char * file, int dmy_fmt);
float mjd(int day, int month, int year, int hour, int minute, int second);
void handle_error(int status);

int main(int argc, char ** argv)
{
    char * inFile  = NULL;
    char * outFile = NULL;
    char * proj    = NULL;
    int in_dmy = 0;
    int fid, dimid;
    int numid, xid, yid, zid, relid, endid;
    int status;
    size_t i;

    // Handle command-line parameters
    for (i = 0; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
            case ARG_OUT:
                if (argc - i < 2) break;
                outFile = argv[++i];
                break;

            case ARG_PROJ:
                if (argc - i < 2) break;
                proj = argv[++i];
                break;

            case ARG_TIME:
                in_dmy = 1;
                break;
            }
        }
        else if (i > 0) inFile = argv[i];
    }
    if (inFile == NULL) {
        printf("Input file not provided!\n%s\n", USAGE_STRING);
        exit(EXIT_FAILURE);
    }
    if (proj == NULL) {
        printf("Projection not provided!\n%s\n", USAGE_STRING);
        exit(EXIT_FAILURE);
    }
    if (outFile == NULL) {
        int len = strlen(inFile);
        outFile = malloc(len);
        strcpy(outFile, inFile);
        outFile[len - 3] = 'n';
        outFile[len - 2] = 'c';
        outFile[len - 1] = '\0';
    }

    // Get and process the particle seed data
    point * cur = read_file(inFile, in_dmy);
    if (proj != NULL) project(proj, cur);

    // Save the data to a new NetCDF file
    create_outfile(outFile, &fid, &numid, &xid, &yid, &zid, &relid, &endid);
    i = 0;
    while (cur != NULL) {
        status = nc_put_var1(fid, numid, &i, &cur->id);
        handle_error(status);
        status = nc_put_var1(fid, xid, &i, &cur->x);
        handle_error(status);
        status = nc_put_var1(fid, yid, &i, &cur->y);
        handle_error(status);
        status = nc_put_var1(fid, zid, &i, &cur->z);
        handle_error(status);
        status = nc_put_var1(fid, relid, &i, &cur->start);
        handle_error(status);
        status = nc_put_var1(fid, endid, &i, &cur->end);
        handle_error(status);
        cur = cur->next;
        i++;
    }

    status = nc_close(fid);
    handle_error(status);
    exit(EXIT_SUCCESS);
}

// Handle status codes returned from the NetCDF library
void handle_error(int status)
{
    if (status != NC_NOERR) {
        fprintf(stderr, "%s\n", nc_strerror(status));
        exit(EXIT_FAILURE);
    }
}

// Format the provided date/time as a Modified Julian Day
float mjd(int day, int month, int year, int hour, int minute, int second)
{
    int a = (14 - month)/12;
    int y = year + 4800 - a;
    int m = month + 12*a - 3;

    int jdn = day + (153*m + 2)/5 + 365*y + y/4 - y/100 + y/400 - 32045;
    double jd = (double)jdn + (double)(hour - 12)/24.0 + (double)minute/1440.0 + (double)second/86400.0;

    return (float)(jd - 2400000.5);
}

// Read the specified particle seeding .dat file into a linked-list of point objects
point * read_file(const char * file, int dmy_fmt)
{
    FILE * fp = fopen(file, "r");
    point * plist = malloc(sizeof(point));
    if (plist == NULL) exit(EXIT_FAILURE);
    point * cur = plist;
    point * prev = NULL;

    while (1) {
        int count = fscanf(fp, "%d %f %f %f", &cur->id, &cur->x, &cur->y, &cur->z);
        if (count == EOF) {
            free(cur);
            if (prev != NULL) prev->next = NULL;
            break;
        }

        if (dmy_fmt == 0) {
            fscanf(fp, "%f %f", &cur->start, &cur->end);
        } else {
            int d, m, y, h, mm, s;
            fscanf(fp, "%d/%d/%d %d:%d:%d", &d, &m, &y, &h, &mm, &s);
            cur->start = mjd(d, m, y, h, mm, s);
            fscanf(fp, "%d/%d/%d %d:%d:%d", &d, &m, &y, &h, &mm, &s);
            cur->end = mjd(d, m, y, h, mm, s);
        }

        prev = cur;
        cur = malloc(sizeof(point));
        if (cur == NULL) exit(EXIT_FAILURE);
        prev->next = cur;
    }
    return plist;
}

// Create a new NetCDF file, formatted for the particle location/release data
void create_outfile(const char * file,
                           int * pfid,
                           int * pnumid,
                           int * pxid,
                           int * pyid,
                           int * pzid,
                           int * prelid,
                           int * pendid)
{
    int status, dimid;

    // Define the variables used in the NetCDF file
    status = nc_create(file, NC_CLOBBER, pfid);
    handle_error(status);
    status = nc_def_dim(*pfid, "number", NC_UNLIMITED, &dimid);
    handle_error(status);
    status = nc_def_var(*pfid, "number", NC_INT, 1, &dimid, pnumid);
    handle_error(status);
    status = nc_def_var(*pfid, "x", NC_FLOAT, 1, &dimid, pxid);
    handle_error(status);
    status = nc_def_var(*pfid, "y", NC_FLOAT, 1, &dimid, pyid);
    handle_error(status);
    status = nc_def_var(*pfid, "z", NC_FLOAT, 1, &dimid, pzid);
    handle_error(status);
    status = nc_def_var(*pfid, "release", NC_FLOAT, 1, &dimid, prelid);
    handle_error(status);
    status = nc_def_var(*pfid, "end", NC_FLOAT, 1, &dimid, pendid);
    handle_error(status);

    // Add descriptive attributes to the new variables
    status = nc_put_att_text(*pfid, *pnumid, "long_name", 19, "Particle Identifier");
    handle_error(status);
    status = nc_put_att_text(*pfid, *pnumid, "standard_name", 2, "id");
    handle_error(status);
    status = nc_put_att_text(*pfid, *pxid, "long_name", 19, "Domain X-Coordinate");
    handle_error(status);
    status = nc_put_att_text(*pfid, *pxid, "units", 6, "meters");
    handle_error(status);
    status = nc_put_att_text(*pfid, *pyid, "long_name", 19, "Domain Y-Coordinate");
    handle_error(status);
    status = nc_put_att_text(*pfid, *pyid, "units", 6, "meters");
    handle_error(status);
    status = nc_put_att_text(*pfid, *pzid, "long_name", 19, "Domain Z-Coordinate");
    handle_error(status);
    status = nc_put_att_text(*pfid, *pzid, "units", 6, "meters");
    handle_error(status);
    status = nc_put_att_text(*pfid, *pzid, "standard_name", 5, "depth");
    handle_error(status);
    status = nc_put_att_text(*pfid, *pzid, "positive", 4, "down");
    handle_error(status);
    status = nc_put_att_text(*pfid, *prelid, "long_name", 21, "Particle Release Time");
    handle_error(status);
    status = nc_put_att_text(*pfid, *prelid, "units", 36, "days since 1858-11-17 00:00:00 +0:00");
    handle_error(status);
    status = nc_put_att_text(*pfid, *pendid, "long_name", 23, "Particle Track End Time");
    handle_error(status);
    status = nc_put_att_text(*pfid, *pendid, "units", 36, "days since 1858-11-17 00:00:00 +0:00");
    handle_error(status);

    // Put the file in data mode, this will also commit the formatting changes we just made
    status = nc_enddef(*pfid);
    handle_error(status);
}

// Perform a coordinate transformation on the point list, to convert spherical coordinates into
// Cartesian coordinates using the provided projection reference
void project(const char * proj, point * plist)
{
    projPJ pj_lcc = pj_init_plus(proj);
    if (pj_lcc == NULL) {
        printf("Failed to initialize projection\n");
        printf("Error: %s\n", pj_strerrno(*pj_get_errno_ref()));
        exit(EXIT_FAILURE);
    }
    projPJ pj_latlong = pj_latlong_from_proj(pj_lcc);
    if (pj_latlong == NULL) {
        printf("Failed to initialize projection\n");
        printf("Error: %s\n", pj_strerrno(*pj_get_errno_ref()));
        exit(EXIT_FAILURE);
    }

    point * cur = plist;
    while (cur != NULL) {
        double x = (double)cur->x * DEG_TO_RAD;
        double y = (double)cur->y * DEG_TO_RAD;
        int status = pj_transform(pj_latlong, pj_lcc, 1, 1, &x, &y, NULL);
        if (status != 0) {
            printf("Failed to perform coordinate projection\n");
            printf("Error %d: %s\n", status, pj_strerrno(status));
            exit(EXIT_FAILURE);
        }
        cur->x = (float)x;
        cur->y = (float)y;
        cur = cur->next;
    }

    pj_free(pj_latlong);
    pj_free(pj_lcc);
}

