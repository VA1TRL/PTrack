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
 *                                                                            *
 * Author: Tristan.Losier@unb.ca                                              *
 ******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <proj_api.h>
#include <netcdf.h>

#define USAGE_STRING  "Command usage:\n    genseed [-o outfile.nc] [-t] [-p projection] particle_release.dat"
#define PROJ_ATT_NAME "CoordinateProjection"

#define ARG_PROJ  'p'
#define ARG_OUT   'o'
#define ARG_TIME  't'
#define ARG_PFILE 'f'

typedef struct _point {
    struct _point * next;
    float x;
    float y;
    float z;
    float start;
    float end;
    int id;
} point;

typedef struct _ncfile {
    int ncid;
    int numberid;
    int xid;
    int yid;
    int zid;
    int releaseid;
    int endid;
} ncfile;

void project(const char * proj, point * plist);
ncfile * create_outfile(const char * file);
point * read_file(const char * file, int dmy_fmt);
char * read_proj(const char * pfile);
float mjd(int day, int month, int year, int hour, int minute, int second);
void handle_error(int status);
void save_points(ncfile * file, point * pts);
void print_points(point * pts);

int main(int argc, char ** argv)
{
    char * inFile  = NULL;
    char * outFile = NULL;
    char * proj    = NULL;
    char * pfile   = NULL;
    int in_dmy = 0;

    // Handle command-line parameters
    int i;
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

            case ARG_PFILE:
                if (argc - i < 2) break;
                pfile = argv[++i];
                break;
            }
        }
        else if (i > 0) inFile = argv[i];
    }

    if (inFile == NULL) {
        printf("Input file not provided!\n%s\n", USAGE_STRING);
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
    point * pts = read_file(inFile, in_dmy);
    if (pfile != NULL) proj = read_proj(pfile);
    if (proj != NULL) {
        printf("Projection reference: %s\n", proj);
        project(proj, pts);
    }

    // Save the data to a new NetCDF file
    ncfile * file = create_outfile(outFile);
    save_points(file, pts);

    nc_close(file->ncid);
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

/* Read the specified particle seeding .dat file into a linked-list of point objects.
 * The .dat file should contain one particle record per line, with the parameters describing
 * the particle separated by one or more spaces. The order and meaning of the parameters
 * that make up the particle record are detailed in the following table.
 *
 *  Column | Name           | Type  | Units                      | Description
 *  ------ | -------------- | ----- | -------------------------- | -----------
 *  1      | ID             | int   | none                       | Arbitrary identifier
 *  2      | X              | float | meters or degrees_east     | Domain x coordinate
 *  3      | Y              | float | meters or degrees_north    | Domain y coordinate
 *  4      | Z              | float | meters                     | Downward positive particle depth
 *  5      | Release Time   | float | MJD or yyyy-mm-dd hh:mm:ss | Time to release particle into the simulation
 *  6      | Track End Time | float | MJD or yyyy-mm-dd hh:mm:ss | Time to remove particle from the simulation
 */
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
            fscanf(fp, "%d-%d-%d %d:%d:%d", &y, &m, &d, &h, &mm, &s);
            cur->start = mjd(d, m, y, h, mm, s);
            fscanf(fp, "%d-%d-%d %d:%d:%d", &y, &m, &d, &h, &mm, &s);
            cur->end = mjd(d, m, y, h, mm, s);
        }

        prev = cur;
        cur = malloc(sizeof(point));
        if (cur == NULL) exit(EXIT_FAILURE);
        prev->next = cur;
    }
    return plist;
}

// Get a projection reference from the provided NetCDF file
char * read_proj(const char * pfile)
{
    char * proj;
    size_t length;
    int status, ncid;

    status = nc_open(pfile, NC_NOWRITE, &ncid);
    handle_error(status);
    status = nc_inq_attlen(ncid, NC_GLOBAL, PROJ_ATT_NAME, &length);
    handle_error(status);
    proj = malloc(length + 1);
    if (proj == NULL) exit(EXIT_FAILURE);

    status = nc_get_att_text(ncid, NC_GLOBAL, PROJ_ATT_NAME, proj);
    handle_error(status);
    proj[length] = '\0';

    status = nc_close(ncid);
    handle_error(status);
    return proj;
}

// Create a new NetCDF file, formatted for the particle location/release data
ncfile * create_outfile(const char * file)
{
    int status, dimid;
    ncfile * nc = malloc(sizeof(ncfile));
    if (nc == NULL) exit(EXIT_FAILURE);

    // Define the variables used in the NetCDF file
    status = nc_create(file, NC_CLOBBER, &nc->ncid);
    handle_error(status);
    status = nc_def_dim(nc->ncid, "number", NC_UNLIMITED, &dimid);
    handle_error(status);
    status = nc_def_var(nc->ncid, "number", NC_INT, 1, &dimid, &nc->numberid);
    handle_error(status);
    status = nc_def_var(nc->ncid, "x", NC_FLOAT, 1, &dimid, &nc->xid);
    handle_error(status);
    status = nc_def_var(nc->ncid, "y", NC_FLOAT, 1, &dimid, &nc->yid);
    handle_error(status);
    status = nc_def_var(nc->ncid, "z", NC_FLOAT, 1, &dimid, &nc->zid);
    handle_error(status);
    status = nc_def_var(nc->ncid, "release", NC_FLOAT, 1, &dimid, &nc->releaseid);
    handle_error(status);
    status = nc_def_var(nc->ncid, "end", NC_FLOAT, 1, &dimid, &nc->endid);
    handle_error(status);

    // Add descriptive attributes to the new variables
    status = nc_put_att_text(nc->ncid, nc->numberid, "long_name", 19, "Particle Identifier");
    handle_error(status);
    status = nc_put_att_text(nc->ncid, nc->numberid, "standard_name", 2, "id");
    handle_error(status);
    status = nc_put_att_text(nc->ncid, nc->xid, "long_name", 19, "Domain X-Coordinate");
    handle_error(status);
    status = nc_put_att_text(nc->ncid, nc->xid, "units", 6, "meters");
    handle_error(status);
    status = nc_put_att_text(nc->ncid, nc->yid, "long_name", 19, "Domain Y-Coordinate");
    handle_error(status);
    status = nc_put_att_text(nc->ncid, nc->yid, "units", 6, "meters");
    handle_error(status);
    status = nc_put_att_text(nc->ncid, nc->zid, "long_name", 19, "Domain Z-Coordinate");
    handle_error(status);
    status = nc_put_att_text(nc->ncid, nc->zid, "units", 6, "meters");
    handle_error(status);
    status = nc_put_att_text(nc->ncid, nc->zid, "standard_name", 5, "depth");
    handle_error(status);
    status = nc_put_att_text(nc->ncid, nc->zid, "positive", 4, "down");
    handle_error(status);
    status = nc_put_att_text(nc->ncid, nc->releaseid, "long_name", 21, "Particle Release Time");
    handle_error(status);
    status = nc_put_att_text(nc->ncid, nc->releaseid, "units", 36, "days since 1858-11-17 00:00:00 +0:00");
    handle_error(status);
    status = nc_put_att_text(nc->ncid, nc->endid, "long_name", 23, "Particle Track End Time");
    handle_error(status);
    status = nc_put_att_text(nc->ncid, nc->endid, "units", 36, "days since 1858-11-17 00:00:00 +0:00");
    handle_error(status);

    // Put the file in data mode, this will also commit the formatting changes we just made
    status = nc_enddef(nc->ncid);
    handle_error(status);
    return nc;
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

// This function accepts a list of particles, and saves the contents of that
// list to the NetCDF file that the provided ncfile object points to.
void save_points(ncfile * file, point * pts)
{
    int status;
    size_t i = 0;
    point * cur = pts;
    while (cur != NULL) {
        status = nc_put_var1(file->ncid, file->numberid, &i, &cur->id);
        handle_error(status);
        status = nc_put_var1(file->ncid, file->xid, &i, &cur->x);
        handle_error(status);
        status = nc_put_var1(file->ncid, file->yid, &i, &cur->y);
        handle_error(status);
        status = nc_put_var1(file->ncid, file->zid, &i, &cur->z);
        handle_error(status);
        status = nc_put_var1(file->ncid, file->releaseid, &i, &cur->start);
        handle_error(status);
        status = nc_put_var1(file->ncid, file->endid, &i, &cur->end);
        handle_error(status);
        cur = cur->next;
        i++;
    }
    printf("Saved %d points\n", (int)i);
}

