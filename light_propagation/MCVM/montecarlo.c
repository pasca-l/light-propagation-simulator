/*

Monte Carlo Program for voxel model

Calculates intensity and partial path length with given settings and model

    Final output is csv file of intensity, ssp, time-resolved ssp,
    but saves binary file of calculation each percent of progress,
    which is added under "../results" relative to the execution file
    (this directory would not generate automatically unlike other data files)

    Settings in settings.conf, model.conf needs to be in the same order,
    however detector and model layer may be added

    Calculation can be resumed from a point of interruption,
    under condition of settings.conf -> is_newfile == 'n'

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

/* files ================================================================== */
#define STRSIZE 50
typedef char string[STRSIZE];

static FILE *fp_base, *fp_model, *fp_note, *fp_data, *fp_pd, *fp_path, *fp_com;
static string baseconf   = "settings.conf";
static string modelconf  = "model.conf";
static string note       = "../data/temporary_note.txt";
static string datafile   = "../data/binary.data";
static string pdfile     = "../data/binary.pd";
static string pathfile   = "../data/binary.ssp";
static string comfile    = "../data/summary.com";

static FILE *fp_int, *fp_pathl, *fp_tpathl, *fp_ssp, *fp_tssp;
static string intfile    = "../data/intensity.csv";
static string pathlfile  = "../data/layer_path.csv";
static string tpathlfile = "../data/layer_tpath.csv";
static string sspfile    = "../data/ssp.csv";
static string tsspfile   = "../data/tssp_map/tssp";
/* ======================================================================== */

/* global variables to be altered ========================================= */
#define MAX_DET     1        /* number of detectors */
#define MAX_LAYER   5        /* number of layers in model*/

#define SCALE       1.0      /* size of pixel*/

#define MAX_X       30       /* data save x (-MAX_X < MAX_X)*/
#define MAX_Y       30       /* data save y (-MAX_Y < MAX_Y)*/
#define MAX_Z       28       /* data save z (0 < MAX_Z)*/

#define MT          400      /* number of division for TR path length */
#define DT          10       /* time(ps) per division */
#define MT_PD       80       /* number of division for TR photon density*/
#define DT_PD       50       /* time(ps) per division */
#define MT_SSP      16       /* number of division for TR SS profile */
#define DT_SSP      250      /* time(ps) per division */
#define INT_CUT     0.0001   /* cutoff intensity for weight */
#define MAX_SCATTER 100000
/* ======================================================================== */

/* global variables as constants ========================================== */
#define VC      0.3  /* speed of light in vacuum (mm/psec) */

#define PI      3.14159265358979324
#define DEG_RAD (PI / 180.0)
#define RAD_DEG (180.0 / PI)
#define COSZERO (1.0 - 1.0E-12)       /* cosine of about 1e-6 rad */

#define MIN_RAND 1.77635683940025e-15  /* =2**-49 (min of randnumber) */

#define TRUE  1
#define FALSE 0

#define N 25
#define M 7
static unsigned long val[N] = {
    0x95f24dab, 0x0b685215, 0xe76ccae7, 0xaf3ec239, 0x715fad23,
    0x24a590ad, 0x69e4b5ef, 0xbf456141, 0x96bc1b7b, 0xa7bdf825,
    0xc1de75b7, 0x8858a9c9, 0x2da87693, 0xb657f9dd, 0xffdc8a9f,
    0x8121da71, 0x8b823ecb, 0x885d05f5, 0x4e20cd47, 0x5a9ad5d9,
    0x512c0c03, 0xea857ccd, 0x4cc1d30f, 0x8891a8a1, 0xa6b7aadb
};
static unsigned long mag[2] = {0x0, 0x8ebfd028};
static int randnum_counter = 0;
/* ======================================================================== */

/* variables to read from file ============================================ */
static char                   is_newfile;  /* if new file (y/n) */
static unsigned long long int phot_input;  /* total number of input photons */
static double                 g;           /* mean cosine of phase function */
static double                 maxpath;     /* max pathlength for cutoff */

static double source_NA;
static double source_minrad;
static double source_maxrad;
static double source_x;
static double source_y;
static double source_z;
static double source_dx;
static double source_dy;
static double source_dz;

static double detector_NA;
static double detector_minrad;
static double detector_maxrad;
static double detector_x[MAX_DET];
static double detector_y[MAX_DET];
static double detector_z[MAX_DET];
static double detector_dx[MAX_DET];
static double detector_dy[MAX_DET];
static double detector_dz[MAX_DET];

static int    width[MAX_LAYER];
static double mus[MAX_LAYER];
static double mua[MAX_LAYER];
static double mut[MAX_LAYER];
/* ======================================================================== */

/* variables for recording ================================================ */
static unsigned long long int phot_start;
static time_t                 time_start;

static unsigned long long int phot_in;        /* input photons */
static unsigned long long int phot_out;       /* output photons */
static unsigned long long int phot_detected;  /* detected photons */
static unsigned long long int phot_overtime;  /* over-time photons */
static unsigned long long int phot_overpath;  /* over-path photons */
static unsigned long long int phot_overscat;  /* over-scattering photons */
static unsigned long long int time_taken = 0; /* time taken */

static double PD[MAX_Z][MAX_X*2+1];
static double TPD[MT_PD][MAX_Z][MAX_X*2+1];
static double SSP[MAX_DET][MAX_Z][MAX_X*2+1][MAX_Y*2+1];
static double TSSP[MAX_DET][MT_SSP][MAX_Z][MAX_X*2+1][MAX_Y*2+1];
static double PPATH[MAX_DET][MAX_LAYER];
static double TPPATH[MAX_DET][MT][MAX_LAYER];
static double INTENSITY[MAX_DET][MT];
static double intensity_for_ssp[MAX_DET][MT_SSP]; /*???*/

/* temporal recording array per input photon */
static double path[MAX_LAYER];
static double path_for_ssp[MAX_Z][MAX_X*2+1][MAX_Y*2+1];
/* ======================================================================== */

/* variables for precalculation =========================================== */
#define TABLEN  10000 /* size of table */
#define REF_IND 1.4

static double stepsize_table[TABLEN];
static double source_costhe[TABLEN];
static double source_sinthe[TABLEN];
static double source_cospsi[TABLEN];
static double source_sinpsi[TABLEN];
static double costhe[TABLEN];
static double sinthe[TABLEN];
static double cospsi[TABLEN];
static double sinpsi[TABLEN];

static int ref_fg;
static double ref_in[1001];
static double ref_out[1001];
/* ======================================================================== */

/* variables for monte carlo part ========================================= */
static int stop_fg; /* flag to stop calculation */
static int walk_fg; /* flag to stop random walk */

static double xold, yold, zold;
static double xnew, ynew, znew;
static double dx, dy, dz;
static int vox, voy, voz;
static double cross_x, cross_y, cross_z; /* coords at border crossing */

static int layer_old, layer_new;

static double step;

static int refx_fg, refy_fg, refz_fg;
static int inref_fg; /* flag for total internal reflection */

static double weight; /* photon weight */
static double eweight; /* photon intensity at exit */

static double totalpath;
static int scatter_count;
/* ======================================================================== */

/* subfunction declaration ================================================ */
static void LoadSettings(); /* loads values from settings file */
static void InitData(); /* initializes variables if file new */
static void LoadData(); /* loads variables if result partially recorded */
static void SaveData(); /* saves variables to binary data; fast save */
static void Summary();
static void SaveDataAsCsv();
static void SetSeed();
static double GenRand();
static void PreCalculatedTable(); /* creates table of precalculated values */

static void MonteCarlo();
    /* function used within MonteCarlo() */
    static void SourcePosition(); /* photon position at source */
    static void SourceDirection(); /* photon direction at source */
    static void CheckLayer(double z); /* layer which photon is in */
    static void NewStepSize(); /* new step size */
    static void NewDirection(); /* new direction by scattering */
    static void RecordVoxelpath(); /* record path at each voxel */
    static void FixPath(); /* fix path length and position */
    static void RecordExit(); /* record property of exited photon */
    static void CalculateRef(); /* new direction by reflection */
/* ======================================================================== */


int main(void){
    LoadSettings();

    if(is_newfile == 'n'){
        LoadData();
    }
    else{
        InitData();
    }

    PreCalculatedTable();
    SetSeed();

    MonteCarlo();

    SaveDataAsCsv();

    exit(0);
}


void LoadSettings(){
    /* error handling for non existing files or content */
    if((fp_base = fopen(baseconf, "r")) == NULL){
        fprintf(stderr, "%s not found or no content\n", baseconf);
        exit(1);
    }
    if((fp_model = fopen(modelconf, "r")) == NULL){
        fprintf(stderr, "%s not found or no content\n", modelconf);
        exit(1);
    }

    /* values of settings.conf */
    fscanf(fp_base, "%c%*[^\n]%*c", &is_newfile);
    fscanf(fp_base, "%llu%*[^\n]%*c", &phot_input);
    fscanf(fp_base, "%lf%*[^\n]%*c", &g);
    fscanf(fp_base, "%d%*[^\n]%*c", &ref_fg);
    fscanf(fp_base, "%lf%*[^\n]%*c", &maxpath);

    /* values of model.conf */
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_NA);
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_minrad);
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_maxrad);
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_x);
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_y);
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_z);
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_dx);
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_dy);
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_dz);

    fscanf(fp_model, "%lf%*[^\n]%*c", &detector_NA);
    fscanf(fp_model, "%lf%*[^\n]%*c", &detector_minrad);
    fscanf(fp_model, "%lf%*[^\n]%*c", &detector_maxrad);
    for(int i = 0; i < MAX_DET; i++){
        fscanf(fp_model, "%lf%*[^\n]%*c", &detector_x[i]);
        fscanf(fp_model, "%lf%*[^\n]%*c", &detector_y[i]);
        fscanf(fp_model, "%lf%*[^\n]%*c", &detector_z[i]);
        fscanf(fp_model, "%lf%*[^\n]%*c", &detector_dx[i]);
        fscanf(fp_model, "%lf%*[^\n]%*c", &detector_dy[i]);
        fscanf(fp_model, "%lf%*[^\n]%*c", &detector_dz[i]);
    }

    for(int i = 0; i < MAX_LAYER; i++){
        fscanf(fp_model, "%d%*[^\n]%*c", &width[i]);
        fscanf(fp_model, "%lf%*[^\n]%*c", &mus[i]);
        fscanf(fp_model, "%lf%*[^\n]%*c", &mua[i]);
        mut[i] = mus[i] + mua[i];
    }

    fclose(fp_base);
    fclose(fp_model);
}

void InitData(){
    phot_start = 0;

    phot_in = 0;
    phot_out = 0;
    phot_detected = 0;
    phot_overtime = 0;
    phot_overpath = 0;
    phot_overscat = 0;

    memset(PPATH, 0, sizeof(PPATH));
    memset(TPPATH, 0, sizeof(TPPATH));
    memset(INTENSITY, 0, sizeof(INTENSITY));
    memset(intensity_for_ssp, 0, sizeof(intensity_for_ssp));

    memset(PD, 0, sizeof(PD));
    memset(TPD, 0, sizeof(TPD));
    memset(SSP, 0, sizeof(SSP));
    memset(TSSP, 0, sizeof(TSSP));
}

void LoadData(){
    /* error handling for non existing files or content */
    if((fp_note = fopen(note, "r")) == NULL){
        fprintf(stderr, "%s not found or no content\n", note);
        exit(1);
    }
    if((fp_data = fopen(datafile, "rb")) == NULL){
        fprintf(stderr, "%s not found or no content\n", datafile);
        exit(1);
    }
    if((fp_pd = fopen(pdfile, "rb")) == NULL){
        fprintf(stderr, "%s not found or no content\n", pdfile);
        exit(1);
    }
    if((fp_path = fopen(pathfile, "rb")) == NULL){
        fprintf(stderr, "%s not found or no content\n", pathfile);
        exit(1);
    }

    /* values of temporary_note.txt */
    fscanf(fp_note, "%llu%*[^\n]%*c", &phot_in);
    fscanf(fp_note, "%llu%*[^\n]%*c", &phot_out);
    fscanf(fp_note, "%llu%*[^\n]%*c", &phot_detected);
    fscanf(fp_note, "%llu%*[^\n]%*c", &phot_overtime);
    fscanf(fp_note, "%llu%*[^\n]%*c", &phot_overpath);
    fscanf(fp_note, "%llu%*[^\n]%*c", &phot_overscat);
    fscanf(fp_note, "%llu%*[^\n]%*c", &time_taken);
    phot_start = phot_in;

    /* values of binary.data */
    fread(PPATH, sizeof(double), MAX_DET*MAX_LAYER, fp_data);
    fread(TPPATH, sizeof(double), MAX_DET*MT*MAX_LAYER, fp_data);
    fread(INTENSITY, sizeof(double), MAX_DET*MT, fp_data);
    fread(intensity_for_ssp, sizeof(double), MAX_DET*MT_SSP, fp_data);

    /* values of binary.pd */
    fread(PD, sizeof(double), MAX_Z*MAX_X*2+1, fp_pd);
    fread(TPD, sizeof(double), MT_PD*MAX_Z*(MAX_X * 2 + 1), fp_pd);

    /* values of binary.ssp */
    fread(SSP, sizeof(double), MAX_DET*MAX_Z*(MAX_X*2+1)*(MAX_Y*2+1), fp_path);
    fread(TSSP, sizeof(double), MAX_DET*MT_SSP*MAX_Z*(MAX_X*2+1)*(MAX_Y*2+1),
          fp_path);

    fclose(fp_note);
    fclose(fp_data);
    fclose(fp_pd);
    fclose(fp_path);
}

void Summary(){
    fp_com = fopen(comfile, "w");
    if (fp_com == NULL){
        fprintf(stderr, "%s can not open\n", comfile);
        exit(1);
    }

    unsigned long long int sec;
    int month, day, hour, min;

    sec = (double)(time(NULL) - time_start) + time_taken;

    day = sec / (60*60*24);
    sec = sec % (60*60*24);
    hour = sec / (60*60);
    sec = sec % (60*60);
    min = sec / 60;
    sec = sec % 60;

    fprintf(fp_com, "<This file is created by 'montecarlo.c'>\n");
    fprintf(fp_com, "input photons        %12llu\n", phot_in);
    fprintf(fp_com, "output photons       %12llu\n", phot_out);
    fprintf(fp_com, "detect photons       %12llu\n", phot_detected);
    fprintf(fp_com, "over time photons    %12llu\n", phot_overtime);
    fprintf(fp_com, "over path photons    %12llu\n", phot_overpath);
    fprintf(fp_com, "over scatter photons %12llu\n", phot_overscat);
    fprintf(fp_com, "time taken : %dD, %dH, %dM, %lluS\n\n",day,hour,min,sec);

    fprintf(fp_com, "layer width mus mua\n");
    for(int i = 0; i < MAX_LAYER; i++){
        fprintf(fp_com, "%d %d %.4f %.4f\n", i, width[i], mus[i], mua[i]);
    }

    fprintf(fp_com, "\nphase function:  g = %.1lf\n", g);
    fprintf(fp_com, "source NA:       %.1f\n", source_NA);
    fprintf(fp_com, "source radius:   Min = %.1f\tMax = %.1f\n",
            source_minrad, source_maxrad);
    fprintf(fp_com, "source:          (x = %.1f, y = %.1f, z = %.1f)\n",
            source_x, source_y, source_z);
    fprintf(fp_com, "detector NA:     %.1f\n", detector_NA);
    fprintf(fp_com, "detector radius: Min = %.1f\tMax = %.1f\n",
            detector_minrad, detector_maxrad);
    for(int i = 0; i < MAX_DET; i++){
        fprintf(fp_com, "detector[%d]:     (x = %.1f, y = %.1f, z = %.1f)\n",
                i, detector_x[i], detector_y[i], detector_z[i]);
    }

    if(stop_fg == TRUE){
        fprintf(fp_com, "\nCOMPLETED");
    }

    fclose(fp_com);
}

void SaveData(){
    /* error handling for non existing files or content */
    if((fp_note = fopen(note, "w")) == NULL){
        fprintf(stderr, "%s can not open\n", note);
        exit(1);
    }
    if((fp_data = fopen(datafile, "wb")) == NULL){
        fprintf(stderr, "%s can not open\n", datafile);
        exit(1);
    }
    if((fp_pd = fopen(pdfile, "wb")) == NULL){
        fprintf(stderr, "%s can not open\n", pdfile);
        exit(1);
    }
    if((fp_path = fopen(pathfile, "wb")) == NULL){
        fprintf(stderr, "%s can not open\n", pathfile);
        exit(1);
    }

    /* values for temporary_note.txt */
    unsigned long long int sec;
    sec = (double)(time(NULL) - time_start) + time_taken;

    fprintf(fp_note, "%12llu\n", phot_in);
    fprintf(fp_note, "%12llu\n", phot_out);
    fprintf(fp_note, "%12llu\n", phot_detected);
    fprintf(fp_note, "%12llu\n", phot_overtime);
    fprintf(fp_note, "%12llu\n", phot_overpath);
    fprintf(fp_note, "%12llu\n", phot_overscat);
    fprintf(fp_note, "%12llu\n", sec);

    /* values for binary.data */
    fwrite(PPATH, sizeof(double), MAX_DET*MAX_LAYER, fp_data);
    fwrite(TPPATH, sizeof(double), MAX_DET*MT*MAX_LAYER, fp_data);
    fwrite(INTENSITY, sizeof(double), MAX_DET*MT, fp_data);
    fwrite(intensity_for_ssp, sizeof(double), MAX_DET*MT_SSP, fp_data);

    /* values for binary.pd */
    int x, y, z, t;

    for(z = 0; z < MAX_Z; z++){
        for(x = 0; x < MAX_X * 2 + 1; x++){
            PD[z][x] /= SCALE * SCALE * SCALE;
        }
    }
    for(t = 0; t < MT_PD; t++){
        for(z = 0; z < MAX_Z; z++){
            for (x = 0; x < MAX_X * 2 + 1; x++){
                TPD[t][z][x] /= SCALE * SCALE * SCALE;
            }
        }
    }
    fwrite(PD, sizeof(double), MAX_Z*MAX_X*2+1, fp_pd);
    fwrite(TPD, sizeof(double), MT_PD*MAX_Z*(MAX_X*2+1), fp_pd);

    /* values for binary.ssp */
    fwrite(SSP, sizeof(double), MAX_DET*MAX_Z*(MAX_X*2+1)*(MAX_Y*2+1), fp_path);
    fwrite(TSSP, sizeof(double), MAX_DET*MT_SSP*MAX_Z*(MAX_X*2+1)*(MAX_Y*2+1),
           fp_path);

    fclose(fp_note);
    fclose(fp_data);
    fclose(fp_pd);
    fclose(fp_path);
}

void SaveDataAsCsv(){
    /* intensity */
    if((fp_int = fopen(intfile, "w")) == NULL){
        fprintf(stderr, "%s can not open\n", intfile);
        exit(1);
    }

    fprintf(fp_int, "time");
    for(int d = 0; d < MAX_DET; d++){
        fprintf(fp_int, ",[det%d]", d);
    }
    fprintf(fp_int, "\n");
    for(int t = 0; t < MT; t++){
        fprintf(fp_int, "%d", (t + 1) * DT);
        for(int d = 0; d < MAX_DET; d++){
            fprintf(fp_int, ",%.6e", INTENSITY[d][t] / phot_in);
        }
        fprintf(fp_int, "\n");
    }
    fclose(fp_int);

    /* path length per layer */
    /* used for low number of photon input to check trend */
    /*
    if((fp_pathl = fopen(pathlfile, "w")) == NULL){
        fprintf(stderr, "%s can not open\n", pathlfile);
        exit(1);
    }

    fprintf(fp_pathl, "layer,path length\n");
    for(int i = 0; i < MAX_LAYER; i++){
        for(int j = 0; j < MAX_DET; j++){
            fprintf(fp_pathl, "%d,%lf\n", i, PPATH[j][i]);
        }
    }
    fclose(fp_pathl);
    */

    /* time-resolved path length per layer */
    if((fp_tpathl = fopen(tpathlfile, "w")) == NULL){
        fprintf(stderr, "%s can not open\n", tpathlfile);
        exit(1);
    }

    fprintf(fp_tpathl, "time");
    for(int d = 0; d < MAX_DET; d++){
        for(int l = 0; l < MAX_LAYER; l++){
            fprintf(fp_tpathl, ",layer%d [det%d]", l, d);
        }
    }
    fprintf(fp_tpathl, "\n");

    for(int t = 0; t < MT; t++){
        fprintf(fp_tpathl, "%d", (t + 1) * DT);
        for(int d = 0; d < MAX_DET; d++){
            for(int l = 0; l < MAX_LAYER; l++){
                if(INTENSITY[d][t] == 0){
                    fprintf(fp_tpathl, "%lf", 0.0);
                }
                else{
                    fprintf(fp_tpathl, ",%lf",
                            TPPATH[d][t][l] / INTENSITY[d][t]);
                }
            }
        }
        fprintf(fp_tpathl, "\n");
    }
    fclose(fp_tpathl);

    /* partial path length (spatial sensitivity profile) */
    if((fp_ssp = fopen(sspfile, "w")) == NULL){
        fprintf(stderr, "%s can not open\n", sspfile);
        exit(1);
    }

    double int_sum = 0;
    for(int d = 0; d < MAX_DET; d++){
        for(int t = 0; t < MT_SSP; t++){
            int_sum += intensity_for_ssp[d][t];
        }

        fprintf(fp_ssp, "[det%d] %d depths\n", d, MAX_Z);
        for(int z = 0; z < MAX_Z; z++){
            for(int y = 0; y < 2*MAX_Y+1; y++){
                for(int x = 0; x < 2*MAX_X+1; x++){
                    if(int_sum == 0){
                        fprintf(fp_ssp, "%lf", 0.0);
                    }
                    else{
                        fprintf(fp_ssp, "%lf", SSP[d][z][x][y] / int_sum);
                    }
                    if(x != 2*MAX_X){
                        fprintf(fp_ssp, ",");
                    }
                }
                fprintf(fp_ssp, "\n");
            }
        }
    }
    fclose(fp_ssp);

    /* time-resolved partial path length */
    for(int t = 0; t < MT_SSP; t++){
        string filename, count;

        sprintf(count, "%d", t);
        strcpy(filename, tsspfile);
        strcat(filename, "(gate=");
        strcat(filename, count);
        strcat(filename, ").csv");

        if((fp_tssp = fopen(filename, "w")) == NULL){
            fprintf(stderr, "%s can not open\n", tsspfile);
            exit(1);
        }

        fprintf(fp_tssp, "Between %d - %d ps\n", t * DT_SSP, (t + 1) * DT_SSP);
        for(int d = 0; d < MAX_DET; d++){
            fprintf(fp_tssp, "[det%d] %d depths\n", d, MAX_Z);
            for(int z = 0; z < MAX_Z; z++){
                for(int y = 0; y < 2*MAX_Y+1; y++){
                    for(int x = 0; x < 2*MAX_X+1; x++){
                        if(intensity_for_ssp[d][t] == 0){
                            fprintf(fp_tssp, "%lf", 0.0);
                        }
                        else{
                            fprintf(fp_tssp, "%lf",
                                TSSP[d][t][z][x][y] / intensity_for_ssp[d][t]);
                        }
                        if(x != 2*MAX_X){
                            fprintf(fp_tssp, ",");
                        }
                    }
                    fprintf(fp_tssp, "\n");
                }
            }
        }
        fclose(fp_tssp);
    }
}

void SetSeed(){
    unsigned long time_now;
    time((time_t *) &time_now); /* sec since 00:00:00 1/1/1970 */

    for(int j = 0; j < N; j++){
        val[j] -= time_now;
    }

    for(int j = 0; j <= 150; j++){ /* exercise the system routine (???) */
        GenRand();
    }
}

double GenRand(){
    unsigned long y;

    if(randnum_counter == N){
        for(int i = 0; i < N - M; i++){
            val[i] = val[i+M] ^ (val[i] >> 1) ^ mag[val[i] % 2];
        }
        for(int i = N - M; i < N; i++){
            val[i] = val[i+M-N] ^ (val[i] >> 1) ^ mag[val[i] % 2];
        }
        randnum_counter = 0;
    }

    y = val[randnum_counter];
    y ^= (y << 7) & 0x2b5b2500;
    y ^= (y << 15) & 0xdb8b0000;
    y &= 0xffffffff;
    y ^= (y >> 16);

    randnum_counter++;

    return ((double) y / (unsigned long) 0xffffffff);
}

void PreCalculatedTable(){
    /* table for step size */
    for(int i = 1; i < TABLEN; i++){
        stepsize_table[i-1] = -log((double) i / TABLEN) / SCALE;
    }

    /* table for direction at source */
    double theta, psi;
    for(int i = 0; i <= TABLEN; i++){
        theta = asin(source_NA) * i / TABLEN;
        psi = 2 * PI * i / (double)TABLEN;

        source_costhe[i] = cos(theta);
        source_cospsi[i] = cos(psi);
        source_sinthe[i] = sqrt(1 - source_costhe[i] * source_costhe[i]);
        source_sinpsi[i] = sqrt(1 - source_cospsi[i] * source_cospsi[i]);
        if(psi > PI){
            source_sinpsi[i] *= -1;
        }
    }

    /* table for reflection */
    double t1, t2, t3, t4, th1, th2;
    if(ref_fg == 0){
        for(int i = 0; i <= 1000; i++){
            ref_in[i] = 0;
            ref_out[i] = 0;
        }
    }
    else if(ref_fg == 1){
        t1 = (REF_IND - 1) / (REF_IND + 1);

        ref_in[0] = 1;
        ref_out[0] = 1;
        ref_in[1000] = t1 * t1;
        ref_out[1000] = ref_in[1000];

        for(int i = 1; i < 1000; i++){
            /* calculating ref_in values */
            t1 = (double) i / 1000;
            t2 = sqrt(1 - t1 * t1);

            th1 = atan2(t2, t1);
            t2 = t2 / REF_IND;
            th2 = atan2(t2, sqrt(1 - t2 * t2));

            t1 = sin(th1 - th2);
            t2 = sin(th1 + th2);
            t3 = cos(th1 - th2);
            t4 = cos(th1 + th2);

            ref_in[i] = 0.5 * t1 * t1 / (t2 * t2) * (1 + t4 * t4 / (t3 * t3));

            /* calculating ref_out values */
            t1 = (double) i / 1000;
            t1 = sqrt(1 - t1 * t1) * REF_IND;
            if(t1 >= 1){
                ref_out[i] = 1;
            }
            else{
                th2 = atan2(t1, sqrt(1 - t1 * t1));
                t1 = sin(th1 - th2);
                t2 = sin(th1 + th2);
                t3 = cos(th1 - th2);
                t4 = cos(th1 + th2);
                ref_out[i] = 0.5 * t1 * t1 / (t2 * t2) * (1 + t4 * t4 /
                                                                    (t3 * t3));
            }
        }
    }
    else if(ref_fg == 2){
        for(int i = 0; i <= 1000; i++){
            ref_in[i] = 0;
            t1 = (double) i / 1000;
            t1 = sqrt(1 - t1 * t1) * REF_IND;
            if(t1 >= 1){
                ref_out[i] = 1;
            }
            else{
                ref_out[i] = 0;
            }
        }
    }
    else{
        printf("error: value of ref_fg not valid\n");
        exit(1);
    }

    /* table for new direction */
    double c1, c2, c3;
    for(int i = 0; i <= TABLEN; i++){
        theta = (double) i / TABLEN;

        if(g == 0){
            costhe[i] = 1 - 2 * theta;
            sinthe[i] = sqrt(1 - costhe[i] * costhe[i]);
        }
        else{
            c1 = 1 + g * g;
            c2 = (1 - g * g) * (1 - g * g);
            c3 = (1 + g - 2 * g * theta) * (1 + g - 2 * g * theta);

            costhe[i] = (c1 - c2 / c3) / (2 * g);
            if(1 - costhe[i] * costhe[i] <= 0){
                sinthe[i] = 0;
            }
            else{
                sinthe[i] = sqrt(1 - costhe[i] * costhe[i]);
            }
        }

        psi = 2 * PI * i / (double)TABLEN;
        cospsi[i] = cos(psi);
        sinpsi[i] = sin(psi);
    }
}


void SourcePosition(){
    double radius, theta;
    radius = source_minrad + ((source_maxrad - source_minrad) * GenRand());
    theta = 2 * PI * GenRand();
    if(source_dx == 1 || source_dx == -1){
        yold = yold + radius * cos(theta);
        zold = zold + radius * sin(theta);
    }
    else if(source_dy == 1 || source_dy == -1){
        xold = xold + radius * cos(theta);
        zold = zold + radius * sin(theta);
    }
    else if(source_dz == 1 || source_dz == -1){
        xold = xold + radius * cos(theta);
        yold = yold + radius * sin(theta);
    }
}

void SourceDirection(){
    long r1, r2;
    double norm, temp;

    r1 = (long) GenRand() * (TABLEN - 1);
    r2 = (long) GenRand() * (TABLEN - 1);

    /* photon input nearly at angle of zero degrees */
    /* normal launch */
    if(fabs(source_dz) > COSZERO){
        dx = source_sinthe[r1] * source_cospsi[r2];
        dy = source_sinthe[r1] * source_sinpsi[r2];
        dz = source_costhe[r1];
        if(source_dz < 0){
            dz *= -1;
        }
    }
    else{
        temp = sqrt(1 - source_dz * source_dz);
        dx = (source_sinthe[r1] * (source_dx * source_dz * source_cospsi[r2] -
             source_dy * source_sinpsi[r2]) / temp +
             source_dx * source_costhe[r1]);
        dy = (source_sinthe[r1] * (source_dy * source_dz * source_cospsi[r2] +
             source_dx * source_sinpsi[r2]) / temp +
             source_dy * source_costhe[r1]);
        dz = (-1 * source_sinthe[r1] * source_cospsi[r2] * temp +
             source_dz * source_costhe[r1]);
    }

    norm = sqrt(dx * dx + dy * dy + dz * dz);
    dx /= norm;
    dy /= norm;
    dz /= norm;
}

void CheckLayer(double z){
    if(z < 0){
        printf("error: invalid z coordinate, (z = %lf)\n", z);
        exit(0);
    }

    int layer_min = 0;
    for(int i = 0; i < MAX_LAYER; i++){
        if(z >= layer_min && z < (layer_min + width[i])){
            layer_new = i;
            break;
        }
        else if(i == MAX_LAYER - 1){
            layer_new = i;
        }
        else{
            layer_min += width[i];
        }
    }
}

void NewStepSize(){
    long r;

    r = (long) (GenRand() * ((double) TABLEN - 2));

    step = stepsize_table[r] / mut[layer_new];
    if(step == 0){
        printf("step size set to V0:\nrand = %ld\nstep = %lf\n",
                r, stepsize_table[r]);
        step = MIN_RAND;
    }
}

void NewDirection(){
    double dr1, dr2, dr3, temp;
    long r1, r2;

    r1 = (long) (GenRand() * ((double) TABLEN - 1));
    r2 = (long) (GenRand() * ((double) TABLEN - 1));

    if(g == 0){
        dx = sinthe[r1] * cospsi[r2];
        dy = sinthe[r1] * sinpsi[r2];
        dz = costhe[r1];
    }
    else{
        dr1 = sinthe[r1] * cospsi[r2];
        dr2 = sinthe[r1] * sinpsi[r2];
        dr3 = costhe[r1];

        temp = 1 - dz * dz;
        if(temp > 0){
            temp = sqrt(temp);
            dx = (dx * dz * dr1 - dy * dr2) / temp + dx * dr3;
            dy = (dy * dz * dr1 + dx * dr2) / temp + dy * dr3;
            dz = dz * dr3 - temp * dr1;
        }
        /* same direction */
        else if(dz > 0){
            dx = dr1;
            dy = dr2;
            dz = dr3;
        }
        /* opposite direction */
        else{
            dx = -1 * dr1;
            dy = -1 * dr2;
            dz = -1 * dr3;
        }
    }
}

void RecordVoxelpath(){
    int fg = 0;
    double lx, ly, lz;
    double temp, xtemp, ytemp, ztemp, ltemp;

    temp = step;
    xtemp = xold;
    ytemp = yold;
    ztemp = zold;
    ltemp = 0;

    /* voxel position */
    if(xtemp < 0){
        vox = (int) xtemp - 1;
    }
    else{
        vox = (int) xtemp;
    }
    if(ytemp < 0){
        voy = (int) ytemp - 1;
    }
    else{
        voy = (int) ytemp;
    }
    voz = (int) ztemp;

    while(fg == 0){
        CheckLayer(ztemp);
        layer_old = layer_new;

        refx_fg = 0;
        refy_fg = 0;
        refz_fg = 0;

        if(dx > 0){
            lx = (vox + 1 - xtemp) / dx;
        }
        else if(dx < 0){
            lx = (vox - xtemp) / dx;
        }
        else if(dx == 0){
            lx = (double) 1 / 0;
        }

        if(dy > 0){
            ly = (voy + 1 - ytemp) / dy;
        }
        else if(dy < 0){
            ly = (voy - ytemp) / dy;
        }
        else if(dy == 0){
            ly = (double) 1 / 0;
        }

        if(dz > 0){
            lz = (voz + 1 - ztemp) / dz;
        }
        else if(dz < 0){
            lz = (voz - ztemp) / dz;
        }
        else if(dz == 0){
            lz = (double) 1 / 0;
        }

        /* possible combinations of photon moving in voxel */
        if(lx < ly && lx < lz && lx <= temp){
            ltemp += lx;
            path[layer_old] += lx;
            if(-1 * MAX_Y <= voy && voy <= MAX_Y &&
               0 <= voz && voz < MAX_Z &&
               -1 * MAX_X <= vox && vox < MAX_X + 1){
                   path_for_ssp[voz][vox+MAX_X][voy+MAX_Y] += lx;
            }

            temp = temp - lx;
            refx_fg = 1;
            xtemp = xtemp + lx * dx;
            ytemp = ytemp + lx * dy;
            ztemp = ztemp + lx * dz;
            if(dx > 0){
                vox += 1;
            }
            else if(dx < 0){
                vox -= 1;
            }
        }
        else if(ly < lx && ly < lz && ly <= temp){
            ltemp += ly;
            path[layer_old] += ly;
            if(-1 * MAX_Y <= voy && voy <= MAX_Y &&
               0 <= voz && voz < MAX_Z &&
               -1 * MAX_X <= vox && vox < MAX_X + 1){
                   path_for_ssp[voz][vox+MAX_X][voy+MAX_Y] += ly;
            }

            temp = temp - ly;
            refy_fg = 1;
            xtemp = xtemp + ly * dx;
            ytemp = ytemp + ly * dy;
            ztemp = ztemp + ly * dz;
            if(dy > 0){
                voy += 1;
            }
            else if(dy < 0){
                voy -= 1;
            }
        }
        else if(lz < lx && lz < ly && lz <= temp){
            ltemp += lz;
            path[layer_old] += lz;
            if(-1 * MAX_Y <= voy && voy <= MAX_Y &&
               0 <= voz && voz < MAX_Z &&
               -1 * MAX_X <= vox && vox < MAX_X + 1){
                   path_for_ssp[voz][vox+MAX_X][voy+MAX_Y] += lz;
            }

            temp = temp - lz;
            refz_fg = 1;
            xtemp = xtemp + lz * dx;
            ytemp = ytemp + lz * dy;
            ztemp = ztemp + lz * dz;
            if(dz > 0){
                voz += 1;
            }
            else if(dz < 0){
                voz -= 1;
            }
        }
        else if(lx == ly && lx < lz && lx <= temp){
            ltemp += lx;
            path[layer_old] += lx;
            if(-1 * MAX_Y <= voy && voy <= MAX_Y &&
               0 <= voz && voz < MAX_Z &&
               -1 * MAX_X <= vox && vox < MAX_X + 1){
                   path_for_ssp[voz][vox+MAX_X][voy+MAX_Y] += lx;
            }

            temp = temp - lx;
            refx_fg = 1;
            refy_fg = 1;
            xtemp = xtemp + lx * dx;
            ytemp = ytemp + lx * dy;
            ztemp = ztemp + lx * dz;
            if(dx > 0){
                vox += 1;
            }
            else if(dx < 0){
                vox -= 1;
            }
            if(dy > 0){
                voy += 1;
            }
            else if(dy < 0){
                voy -= 1;
            }
        }
        else if(lx < ly && lx == lz && lx <= temp){
            ltemp += lx;
            path[layer_old] += lx;
            if(-1 * MAX_Y <= voy && voy <= MAX_Y &&
               0 <= voz && voz < MAX_Z &&
               -1 * MAX_X <= vox && vox < MAX_X + 1){
                   path_for_ssp[voz][vox+MAX_X][voy+MAX_Y] += lx;
            }

            temp = temp - lx;
            refx_fg = 1;
            refz_fg = 1;
            xtemp = xtemp + lx * dx;
            ytemp = ytemp + lx * dy;
            ztemp = ztemp + lx * dz;
            if(dx > 0){
                vox += 1;
            }
            else if(dx < 0){
                vox -= 1;
            }
            if(dz > 0){
                voz += 1;
            }
            else if(dz < 0){
                voz -= 1;
            }
        }
        else if(ly < lx && ly == lz && ly <= temp){
            ltemp += ly;
            path[layer_old] += ly;
            if(-1 * MAX_Y <= voy && voy <= MAX_Y &&
               0 <= voz && voz < MAX_Z &&
               -1 * MAX_X <= vox && vox < MAX_X + 1){
                   path_for_ssp[voz][vox+MAX_X][voy+MAX_Y] += ly;
            }

            temp = temp - ly;
            refy_fg = 1;
            refz_fg = 1;
            xtemp = xtemp + ly * dx;
            ytemp = ytemp + ly * dy;
            ztemp = ztemp + ly * dz;
            if(dy > 0){
                voy += 1;
            }
            else if(dy < 0){
                voy -= 1;
            }
            if(dz > 0){
                voz += 1;
            }
            else if(dz < 0){
                voz -= 1;
            }
        }
        else if(lx == ly && ly == lz && lx == lz && lx <= temp){
            ltemp += lx;
            path[layer_old] += lx;
            if(-1 * MAX_Y <= voy && voy <= MAX_Y &&
               0 <= voz && voz < MAX_Z &&
               -1 * MAX_X <= vox && vox < MAX_X + 1){
                   path_for_ssp[voz][vox+MAX_X][voy+MAX_Y] += lx;
            }

            temp = temp - lx;
            refx_fg = 1;
            refy_fg = 1;
            refz_fg = 1;
            xtemp = xtemp + lx * dx;
            ytemp = ytemp + lx * dy;
            ztemp = ztemp + lx * dz;
            if(dx > 0){
                vox += 1;
            }
            else if(dx < 0){
                vox -= 1;
            }
            if(dy > 0){
                voy += 1;
            }
            else if(dy < 0){
                voy -= 1;
            }
            if(dz > 0){
                voz += 1;
            }
            else if(dz < 0){
                voz -= 1;
            }
        }
        else if(temp < lx && temp < ly && temp < lz){
            ltemp += temp;
            path[layer_old] += temp;
            if(-1 * MAX_Y <= voy && voy <= MAX_Y &&
               0 <= voz && voz < MAX_Z &&
               -1 * MAX_X <= vox && vox < MAX_X + 1){
                   path_for_ssp[voz][vox+MAX_X][voy+MAX_Y] += temp;
            }

            temp = 0;
            fg = 1;
        }
        else{
            printf("error: unexpected condition given in RecordVoxelpath\n");
            printf("(lx,ly,lz,temp) = (%lf,%lf,%lf,%lf)\n", lx, ly, lz, temp);
            exit(1);
        }

        if(voz < 0){
            ltemp += temp;
            fg = 1;
            ztemp = 0;
            cross_x = xtemp;
            cross_y = ytemp;
            cross_z = ztemp;
        }
        if(ztemp < 0){
            printf("error: invalid z coord\nz = %lf\nvoz = %d\n", ztemp, voz);
            exit(1);
        }

        CheckLayer(ztemp);
        if(layer_new != layer_old){
            temp = temp * mut[layer_old] / mut[layer_new];
        }
    }
    step = ltemp;
}

void FixPath(){
    double temp;

    temp = sqrt((xnew - cross_x) * (xnew - cross_x) +
                (ynew - cross_y) * (ynew - cross_y) +
                (znew - cross_z) * (znew - cross_z));
    totalpath -= temp;
    step = temp;

    if(step == 0){
        step = -1 * log((double) MIN_RAND) / mut[0] / SCALE;
    }

    xnew = cross_x;
    ynew = cross_y;
    znew = cross_z;
}

void RecordExit(){
    int x, y, z;
    long t, ssp_t;
    double theta, d1, d2, temp;

    t = (long) fabs((totalpath * REF_IND * SCALE / VC) / DT);
    ssp_t = (long) fabs((totalpath * REF_IND * SCALE / VC) / DT_SSP);

    for(int i = 0; i < MAX_DET; i++){
        theta = acos(detector_dx[i] * dx +
                     detector_dy[i] * dy +
                     detector_dz[i] * dz);

        if(detector_dz[i] == 1 || detector_dz[i] == -1){
            d1 = (detector_x[i] - xnew) * (detector_x[i] - xnew);
            d2 = (detector_y[i] - ynew) * (detector_y[i] - ynew);
        }

        if(sqrt(d1 + d2) >= (double) detector_minrad / SCALE &&
           sqrt(d1 + d2) <= (double) detector_maxrad / SCALE &&
           theta <= asin(detector_NA)){
               phot_detected++;

               for(int j = 0; j < MAX_LAYER; j++){
                   PPATH[i][j] += eweight * path[j] * SCALE;
                   if(t >= 0 && t < MT){
                       TPPATH[i][t][j] += eweight * path[j] * SCALE;
                   }
               }
               for(z = 0; z < MAX_Z; z++){
                   for(x = 0; x < 2 * MAX_X; x++){
                       for(y = 0; y <= 2 * MAX_Y; y++){
                           temp = eweight * path_for_ssp[z][x][y] * SCALE;
                           SSP[i][z][x][y] += temp;
                       }
                   }
               }

               if(ssp_t >= 0 && ssp_t < MT_SSP){
                   intensity_for_ssp[i][ssp_t] += eweight;
                   for(z = 0; z < MAX_Z; z++){
                       for(x = 0; x < 2 * MAX_X; x++){
                           for(y = 0; y <= 2 * MAX_Y; y++){
                               temp = eweight * path_for_ssp[z][x][y] * SCALE;
                               TSSP[i][ssp_t][z][x][y] += temp;
                           }
                       }
                   }
               }

               if(t >= 0 && t < MT){
                   INTENSITY[i][t] += eweight;
               }
        }
    }
}

void CalculateRef(){
    double n1, n2, n3, ref, temp;

    if(refx_fg == 1 && refy_fg == 0 && refz_fg == 0){
        n1 = -1 * dx / fabs(dx);
        n2 = 0;
        n3 = 0;
    }
    else if(refx_fg == 0 && refy_fg == 1 && refz_fg == 0){
        n1 = 0;
        n2 = -1 * dy / fabs(dy);
        n3 = 0;
    }
    else if(refx_fg == 0 && refy_fg == 0 && refz_fg == 1){
        n1 = 0;
        n2 = 0;
        n3 = -1 * dy / fabs(dz);
    }
    else{
        n1 = -1 * dx;
        n2 = -1 * dy;
        n3 = -1 * dz;
    }

    temp = (-1 * dx * n1 + -1 * dy * n2 + -1 * dz * n3);
    dx = 2 * temp * n1 + dx;
    dy = 2 * temp * n2 + dy;
    dz = 2 * temp * n3 + dz;
    ref = ref_out[(int) (fabs(temp * 1000) + 0.5)];

    if(ref == 1 || 1 - (1 - temp * temp) * REF_IND * REF_IND < 0){
        inref_fg = 1;
    }
    else{
        eweight = weight * (1 - ref);
        weight = weight * ref;
        inref_fg = 0;
    }
}

void MonteCarlo(){
    long t, time_max;
    double pd;

    time_start = time(NULL);
    time_max = MT * DT;
    stop_fg = FALSE;

    while(phot_in < phot_input && stop_fg == FALSE){
        phot_in ++;
        // path and path_for_ssp might be only for temporary recording
        // which needs to for SSP to be determined
        memset(path, 0, sizeof(path));
        memset(path_for_ssp, 0, sizeof(path_for_ssp));
        totalpath = 0;

        xold = source_x;
        yold = source_y;
        zold = source_z;
        SourcePosition();
        SourceDirection();

        if(zold < 0){
            printf("error: invalid source position\n");
            exit(1);
        }
        CheckLayer(zold);
        NewStepSize();

        weight = 1;
        eweight = 0;

        inref_fg = 0;
        walk_fg = TRUE;
        scatter_count = 0;
        while(walk_fg == TRUE){
            scatter_count++;
            RecordVoxelpath();
            totalpath += step;
            xnew = xold + step * dx;
            ynew = yold + step * dy;
            znew = zold + step * dz;
            t = (long)(totalpath * REF_IND * SCALE / VC);

            /* when photon exits */
            if(znew < 0){
                FixPath();
                if (ref_fg == 1){
                    CalculateRef();
                }
                else{
                    eweight = weight;
                    weight = 0;
                    dx = -1 * dx;
                    dy = -1 * dy;
                    dz = -1 * dz;
                }

                if(inref_fg == 0){
                    phot_out++;
                    RecordExit();
                    if (weight < INT_CUT){
                        walk_fg = FALSE;
                    }
                }
            }
            else if(t > time_max){
                phot_overtime++;
                walk_fg = FALSE;
            }
            else if(totalpath * SCALE >= maxpath){
                phot_overpath++;
                walk_fg = FALSE;
            }
            else if(scatter_count >= MAX_SCATTER-1){
                phot_overscat++;
                walk_fg = FALSE;
            }
            else{
                xold = xnew;
                yold = ynew;
                zold = znew;

                if (voy == (int)source_y &&
                    -1 * MAX_X <= vox && vox < MAX_X + 1 &&
                    voz < MAX_Z){
                        pd = weight / mut[layer_new] *
                             (mut[layer_new] - mus[layer_new]);
                        PD[voz][vox + MAX_X] += pd;
                }

                if (voy == -1 * MAX_Y &&
                    -1 * MAX_X <= vox && vox < MAX_X + 1 &&
                    voz < MAX_Z){
                        if ((int)(t / DT_PD) < MT_PD){
                            pd = weight / mut[layer_new] *
                                 (mut[layer_new] - mus[layer_new]);
                            TPD[(int)(t / DT_PD)][voz][vox + MAX_X] += pd;
                        }
                }

                weight = weight * mus[layer_new] / mut[layer_new];

                NewStepSize();
                NewDirection();
            }
        }
        if(phot_in >= phot_input){
            stop_fg = TRUE;
        }
        if(phot_in % (phot_input / 1000) == 0){
            Summary();
        }
        if(phot_in % (phot_input / 100) == 0){
            SaveData();

            /* show progress */
            printf("%lld%% done\n", phot_in / (phot_input / 100));
        }
    }
}
