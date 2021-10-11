// https://omlc.org/classroom/ece532/class4/index.html

// minimal montecarlo
// https://omlc.org/classroom/ece532/class4/trmc/trmc.c


/*

written by S.Yamadate[2021]
Monte Carlo Program for voxel model
Put in effort on comments for readability

adapted from mc_fnirs_slab.c written by K.Takai[2014]

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

/* files =================================================================== */
#define STRSIZE 50
typedef char string[STRSIZE];

static FILE   *fp_base, *fp_model, *fp_note, *fp_data, *fp_pd, *fp_ssp, *fp_com;
static string baseconf = "settings.conf";
static string modelconf = "model.conf";
static string note = "temporary_note.txt";
static string datafile = "binary.data";
static string pdfile = "binary.pd";
static string sspfile = "binary.ssp";
static string comfile = "result.com";
/* ========================================================================= */

/* global variables to be altered ========================================== */
#define SCALE       1.0      /* size of pixel*/

#define MAX_X       60       /* data save x (-MAX_X < MAX_X)*/
#define MAX_Y       2        /* data save x (-MAX_X < MAX_X)*/
#define MAX_Z       50       /* data save z (0 < MAX_Z)*/

#define PD_MAX_LINE 5
#define PD_MIN_LINE 1        /*y�������̋��߂�PD�̍��W�͈́i-80�`80�j�i�͈͂�10�ȓ��j*/

#define MAX_DET     1
#define MAX_LAYER   5

#define MT          400      /* number of division for TPSF */
#define DT          10
#define MT_PD       80       /* number of division for TR photon density*/
#define DT_PD       50
#define MT_SSP      16       /* number of division for TR spatial */
    /*sencitivity profile*/
#define DT_SSP      250
#define INT_CUT     0.0001   /* cutoff intensity for cweight */
#define MAX_SCATTER 100000
/* ========================================================================= */

/* global variables as constants =========================================== */
#define VC          0.3  /* speed of light in vacuum (mm/psec) */

#define PI          3.14159265358979324
#define DEG_RAD     (PI / 180.0)
#define RAD_DEG     (180.0 / PI)
#define COSZERO     (1.0 - 1.0E-12)       /* cosine of about 1e-6 rad */

#define V0          1.77635683940025e-15  /* =2**-49 (min of randnumber) */

#define TRUE        1
#define FALSE       0

// extern time_t time(time_t *);
#define N           25
#define M           7
static unsigned long x[N] = {
    0x95f24dab, 0x0b685215, 0xe76ccae7, 0xaf3ec239, 0x715fad23,
    0x24a590ad, 0x69e4b5ef, 0xbf456141, 0x96bc1b7b, 0xa7bdf825,
    0xc1de75b7, 0x8858a9c9, 0x2da87693, 0xb657f9dd, 0xffdc8a9f,
    0x8121da71, 0x8b823ecb, 0x885d05f5, 0x4e20cd47, 0x5a9ad5d9,
    0x512c0c03, 0xea857ccd, 0x4cc1d30f, 0x8891a8a1, 0xa6b7aadb
};
static unsigned long mag[2] = {0x0, 0x8ebfd028};
static int k = 0;
/* ========================================================================= */

/* variables to read from file ============================================= */
static char                   is_newfile;  /* if new file (y/n) */
static unsigned long long int phot_input;  /* total number of input photons */
static double                 g;           /* mean cosine of phase function */
static double                 maxpath;     /* max pathlength for cutoff */

static double source_NA;
static double source_minradius;
static double source_maxradius;
static double source_x;
static double source_y;
static double source_z;
static double source_dx;
static double source_dy;
static double source_dz;

static double detector_NA;
static double detector_minradius;
static double detector_maxradius;
static double detector_x[MAX_DET];
static double detector_y[MAX_DET];
static double detector_z[MAX_DET];
static double detector_dx[MAX_DET];
static double detector_dy[MAX_DET];
static double detector_dz[MAX_DET];

static int    thickness[MAX_LAYER];
static double mus[MAX_LAYER];
static double mua[MAX_LAYER];
static double mut[MAX_LAYER];
/* ========================================================================= */

/* variables for recording ================================================= */
static unsigned long long int phot_start, phot_end;
static time_t                 time_start, time_end;

static unsigned long long int phot_in;        /* input photons */
static unsigned long long int phot_out;       /* output photons */
static unsigned long long int phot_detected;  /* detected photons */
static unsigned long long int phot_overtime;  /* over-time photons */
static unsigned long long int phot_overpath;  /* over-path photons */
static unsigned long long int phot_overscat;  /* over-scattering photons */

static double PD[MAX_Z][MAX_X*2+1];
static double TPD[MT_PD][MAX_Z][MAX_X * 2 + 1];
static double STPD[MT_PD];
static double SSP[MAX_DET][MAX_Z][MAX_X*2+1][MAX_Y*2+1];
static double TSSP[MAX_DET][MT_SSP][MAX_Z][MAX_X*2+1][MAX_Y*2+1];

static double path[MAX_LAYER];
static double path_y[MAX_LAYER];
static double path_for_ssp[MAX_Z][MAX_X*2+1][MAX_Y * 2 + 1];
static double ppath[MAX_DET][MAX_LAYER];
static double tppath[MAX_DET][MT][MAX_LAYER];
static double tppath_y[MAX_DET][MT][MAX_LAYER];
static double intensity[MAX_DET][MT];
static double intensity_for_ssp[MAX_DET][MT_SSP];
/* ========================================================================= */

/* variables for precalculation ============================================ */
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
/* ========================================================================= */


/* variables for monte carlo part ========================================== */
static int	stop_fg;     /* flag to stop calculation */
/* ========================================================================= */



/*unsorted variables*/
/* ========================================================================= */
#define ref_ind		1.4

static double	allpath;
static double	NA_Cos[TABLEN];	/* lookup table for direction cosins at source*/
static double	NA_Sin[TABLEN];	/* lookup table for direction cosins at source*/
static double	intensity[MAX_DET][MT];
static double	intensity_for_ssp[MAX_DET][MT_SSP];

static int	inr;		/* flag for total internal reflection */
static int	fg;			/* flag to stop random walk */
static int	source_fg;
static int	refx_fg;
static int	refy_fg;
static int	refz_fg;


static double 	xold, yold, zold;	/* old x,y,z coordinates */
static double	xnew, ynew, znew;	/* new x,y,z coordinates */
static int		vox,voy,voz;

static int		source_time;       /*[20180209�ǈ�]�Ǝˎ��Ԃ̈ʒu*/

static double	cross_x, cross_y, cross_z;
static int		layer_old,layer_new;
static int		n_scatter;
static double	Len;	/* free path lenght */
static double	Lsub;
static double	dx, dy, dz;	/* direction cosines */
static double	sweight;		/* photon intensity for exit (mua=0) */
static double	cweight;		/* photon intensity for continuation */

static string angfile;
/* ========================================================================= */
/* ========================================================================= */


/* subfunction declaration ================================================= */
static void LoadSettings(); /* loads values from settings file */
static void InitData(); /* initializes variables if file new */
static void LoadData(); /* loads variables if result partially recorded */
static void SaveData(); /* saves variables to binary data */
static void SetSeed();
static double GenRand();
static void PreCalculatedTable(); /* creates table of precalculated values */

static void Output_ref();
static void Output_no_ref();
static void Source_pos();
static void Source_dir();
static void Source_time();
static void New_len();
static void New_dir();
static void Check_layer(double z);
static void Restore_voxelpath();
static void Restore_data();
static void Store_data();
static void Fnstop();
static void Main_mc();
/* ========================================================================= */

/* dummy variables ========================================================= */
long i, j;
double temp;
static char ch = 'n';
/* ========================================================================= */


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

    time_start = time(NULL);
    Main_mc();
    Fnstop();
    SaveData();

    exit(0);
}


void LoadSettings(){
    /* error handling for non existing files or content */
    if((fp_base = fopen(baseconf, "r")) == NULL){
        fprintf(stderr, "%s not found or no content\n", baseconf);
        exit(1);
    }
    if((fp_model = fopen(baseconf, "r")) == NULL){
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
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_minradius);
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_maxradius);
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_x);
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_y);
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_z);
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_dx);
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_dy);
    fscanf(fp_model, "%lf%*[^\n]%*c", &source_dz);

    fscanf(fp_model, "%lf%*[^\n]%*c", &detector_NA);
    fscanf(fp_model, "%lf%*[^\n]%*c", &detector_minradius);
    fscanf(fp_model, "%lf%*[^\n]%*c", &detector_maxradius);
    for(i = 0; i < MAX_DET; i++){
        fscanf(fp_model, "%lf%*[^\n]%*c", &detector_x[i]);
        fscanf(fp_model, "%lf%*[^\n]%*c", &detector_y[i]);
        fscanf(fp_model, "%lf%*[^\n]%*c", &detector_z[i]);
        fscanf(fp_model, "%lf%*[^\n]%*c", &detector_dx[i]);
        fscanf(fp_model, "%lf%*[^\n]%*c", &detector_dy[i]);
        fscanf(fp_model, "%lf%*[^\n]%*c", &detector_dz[i]);
    }

    for(i = 0; i < MAX_LAYER; i++){
        fscanf(fp_model, "%d%*[^\n]%*c", &thickness[i]);
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

    memset(ppath, 0, sizeof(ppath));
    memset(tppath, 0, sizeof(tppath));
    memset(tppath_y, 0, sizeof(tppath_y));
    memset(intensity, 0, sizeof(intensity));
    memset(intensity_for_ssp, 0, sizeof(intensity_for_ssp));

    memset(PD, 0, sizeof(PD));
    memset(TPD, 0, sizeof(TPD));
    memset(STPD, 0, sizeof(STPD));
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
    if((fp_ssp = fopen(sspfile, "rb")) == NULL){
        fprintf(stderr, "%s not found or no content\n", sspfile);
        exit(1);
    }

    /* values of temporary_note.txt */
    fscanf(fp_note, "%llu%*[^\n]%*c", &phot_in);
    fscanf(fp_note, "%llu%*[^\n]%*c", &phot_out);
    fscanf(fp_note, "%llu%*[^\n]%*c", &phot_detected);
    fscanf(fp_note, "%llu%*[^\n]%*c", &phot_overtime);
    fscanf(fp_note, "%llu%*[^\n]%*c", &phot_overpath);
    fscanf(fp_note, "%llu%*[^\n]%*c", &phot_overscat);
    phot_start = phot_in;

    /* values of binary.data */
    fread(ppath, sizeof(double), MAX_DET*MAX_LAYER, fp_data);
    fread(tppath, sizeof(double), MAX_DET*MT*MAX_LAYER, fp_data);
    fread(tppath_y, sizeof(double), MAX_DET*MT*MAX_LAYER, fp_data);
    fread(intensity, sizeof(double), MAX_DET*MT, fp_data);
    fread(intensity_for_ssp, sizeof(double), MAX_DET*MT_SSP, fp_data);

    /* values of binary.pd */
    fread(PD, sizeof(double), MAX_Z*MAX_X*2+1, fp_pd);
    fread(TPD, sizeof(double), MT_PD*MAX_Z*(MAX_X * 2 + 1), fp_pd);
    fread(STPD, sizeof(double), 4000+1, fp_pd);

    /* values of binary.ssp */
    fread(SSP, sizeof(double), MAX_DET*MAX_Z*(MAX_X*2+1)*(MAX_Y*2+1), fp_ssp);
    fread(TSSP, sizeof(double), MAX_DET*MT_SSP*MAX_Z*(MAX_X*2+1)*(MAX_Y*2+1),
          fp_ssp);

    fclose(fp_note);
    fclose(fp_data);
    fclose(fp_pd);
    fclose(fp_ssp);
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
    if((fp_ssp = fopen(sspfile, "wb")) == NULL){
        fprintf(stderr, "%s can not open\n", sspfile);
        exit(1);
    }

    /* values for temporary_note.txt */
    fprintf(fp_note, "%12llu\n", phot_in);
    fprintf(fp_note, "%12llu\n", phot_out);
    fprintf(fp_note, "%12llu\n", phot_detected);
    fprintf(fp_note, "%12llu\n", phot_overtime);
    fprintf(fp_note, "%12llu\n", phot_overpath);
    fprintf(fp_note, "%12llu\n", phot_overscat);

    /* values for binary.data */
    fwrite(ppath, sizeof(double), MAX_DET*MAX_LAYER, fp_data);
    fwrite(tppath, sizeof(double), MAX_DET*MT*MAX_LAYER, fp_data);
    fwrite(tppath_y, sizeof(double), MAX_DET*MT*MAX_LAYER, fp_data);
    fwrite(intensity, sizeof(double), MAX_DET*MT, fp_data);
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
    for(t = 0; t < 4000; t++){
        STPD[t] /= SCALE * SCALE * SCALE;
    }
    fwrite(PD, sizeof(double), MAX_Z*MAX_X*2+1,fp_pd);
    fwrite(TPD, sizeof(double), MT_PD*MAX_Z*(MAX_X * 2 + 1), fp_pd);
    fwrite(STPD, sizeof(double), 4000+1, fp_pd);

    /* values for binary.ssp */
    fwrite(SSP, sizeof(double), MAX_DET*MAX_Z*(MAX_X*2+1)*(MAX_Y*2+1), fp_ssp);
    fwrite(TSSP, sizeof(double), MAX_DET*MT_SSP*MAX_Z*(MAX_X*2+1)*(MAX_Y*2+1),
           fp_ssp);

    fclose(fp_note);
    fclose(fp_data);
    fclose(fp_pd);
    fclose(fp_ssp);
}

void SetSeed(){
    unsigned long time_now;
    time((time_t *) &time_now); /* sec since 00:00:00 1/1/1970 */

    for(j = 0; j < N; j++){
        x[j] -= time_now;
    }

    for(j = 0; j <= 150; j++){ /* exercise the system routine (???) */
        GenRand();
    }
}

double GenRand(){
    unsigned long y;

    if(k == N){
        for(i = 0; i < N - M; i++){
            x[i] = x[i+M] ^ (x[i] >> 1) ^ mag[x[i] % 2];
        }
        for(i = N - M; i < N; i++){
            x[i] = x[i+M-N] ^ (x[i] >> 1) ^ mag[x[i] % 2];
        }
        k = 0;
    }

    y = x[k];
    y ^= (y << 7) & 0x2b5b2500;
    y ^= (y << 15) & 0xdb8b0000;
    y &= 0xffffffff;
    y ^= (y >> 16);

    k++;

    return ((double) y / (unsigned long) 0xffffffff);
}

void PreCalculatedTable(){
    /* table for step size */
    for(i = 1; i < TABLEN; i++){
        stepsize_table[i-1] = -log((double) i / TABLEN) / SCALE;
    }

    /* table for direction at source */
    double theta, psi;
    for(i = 0; i <= TABLEN; i++){
        theta = asin(source_NA) * i / (double)TABLEN;
        psi = 2 * PI * i / (double)TABLEN;

        source_costhe[i] = cos(theta);
        source_cospsi[i] = cos(psi);
        source_sinthe[i] = sqrt(1 - source_costhe[i] * source_costhe[i]);
        source_sinpsi[i] = sqrt(1 - source_cospsi[i] * source_cospsi[i]);
        if(psi > PI){
            source_sinpsi[i] = -1 * source_sinpsi[i];
        }
    }

    /* table for new direction */
    double c1, c2, c3;
    for(i = 0; i <= TABLEN; i++){
        j = (double) i / TABLEN;

        if(g == 0){
            costhe[i] = 1 - 2 * j;
            sinthe[i] = sqrt(1 - costhe[i] * costhe[i]);
        }
        else{
            c1 = 1 + g * g;
            c2 = (1 - g * g) * (1 - g * g);
            c3 = (1 + g - 2 * g * j) * (1 + g - 2 * g * j);

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

    /* table for reflection */
    double t1, t2, t3, t4, th1, th2;
    if(ref_fg == 0){
        for(i = 0; i <= 1000; i++){
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

        for(i = 1; i < 1000; i++){
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
        for(i = 0; i <= 1000; i++){
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
}


static void Main_mc(){
	long time,time_max;
	double pd;

	time_max = MT * DT;

	stop_fg = FALSE;

	Fnstop();
	do{
		phot_in ++;
		memset(path, 0, sizeof(path));
		memset(path_y, 0, sizeof(path_y));
		memset(path_for_ssp, 0, sizeof(path_for_ssp));
		allpath=0.0;
		source_fg = 0;

		/*�� Source�̐ݒ肱������ ��*/
		xold = source_x;
		yold = source_y;
		zold = source_z;
		Source_pos();	//xold,yold,zold�̌���
		Source_dir();	//du�̌���

		Source_time(); //�Ǝˎ��Ԃ̌���

		if(zold < 0){
			printf("error in source\n");
			printf("new position(x,y,z)=(%f,%f,%f)\n",xold,yold,zold);
			exit(1);
		}
		Check_layer(zold);
		New_len();	//Len�̌���

		/*�� Source�̐ݒ肱���܂� ��*/

		cweight = 1.0;
		sweight = 0.0;
		inr=0;
		fg = FALSE;
		n_scatter=0;
		do{
			n_scatter += 1;
			Restore_voxelpath();
			allpath += Len;
			xnew = xold+Len*dx;
			ynew = yold+Len*dy;
			znew = zold+Len*dz;
			time = (long)(allpath*ref_ind*SCALE/VC + source_time);
			if (znew < 0){	/* photon exits */
				Restore_data();		/* correct coordinates */
				if (ref_fg == 1){
					Output_ref();
				}
				else{
					Output_no_ref();
				}
				if(inr == 0){
					phot_out++;
					Store_data();
					if (cweight<INT_CUT){
						fg = TRUE;
					}
				}


				xold = xnew;
				yold = ynew;
				zold = znew;
			}
			else if (time > time_max){
				phot_overtime++;
				fg = TRUE;
			}
			else if (allpath*SCALE>=maxpath){	/* lost photon (cut off) */
				phot_overpath++;
				fg = TRUE;
			}
			else if (n_scatter >= MAX_SCATTER-1){
				phot_overscat++;
				fg = TRUE;
			}
			else{

				xold = xnew;
				yold = ynew;
				zold = znew;

				if ((voy == (int)source_y) && (vox >= -MAX_X) && (vox < MAX_X + 1) && (voz < MAX_Z)){
					pd = cweight / mut[layer_new] * (mut[layer_new] - mus[layer_new]);
					PD[voz][vox + MAX_X] += pd;
				}


			if ((PD_MIN_LINE == voy) && (vox >= -MAX_X) && (vox < MAX_X + 1) && (voz < MAX_Z)){
					if ((int)(time / DT_PD) < MT_PD){
						pd = cweight / mut[layer_new] * (mut[layer_new] - mus[layer_new]);
						TPD[(int)(time / DT_PD)][voz][vox + MAX_X] += pd;
					}
				}


				cweight = cweight*mus[layer_new]/mut[layer_new];

				New_len();
				New_dir();
			}
		} while (fg != TRUE);
		if (phot_in % 100 == 0){
			Fnstop();
		}
		if (phot_in % 10000000 == 0){
			SaveData();
		}
        if (phot_in < phot_input){
            ch = 'y';
        }
	} while ((phot_in < phot_input) && (stop_fg!=TRUE));

}

static void Output_ref(){
	/* calculate reflection and refraction in going from tissue to air. */
	double refl, temp, n_ref1, n_ref2, n_ref3, trans_ratio, theta, duu, ds1, ds2, ds3;

	duu = 0.0;
	refx_fg = 0;
	refy_fg = 0;
	refz_fg = 1;

	if(refx_fg == 1 && refy_fg == 0 && refz_fg == 0){
	n_ref1 = -1.0*dx/fabs(dx);
	n_ref2 = 0.0;
	n_ref3 = 0.0;
	//printf("n_ref1 = %lg \n",n_ref1);
	}
	else if(refx_fg == 0 && refy_fg == 1 && refz_fg == 0){
		n_ref1 = 0.0;
		n_ref2 = -1.0*dy/fabs(dy);
		n_ref3 = 0.0;
		//printf("n_ref2 = %lg \n",n_ref2);
	}
	else if(refx_fg == 0 && refy_fg ==0 && refz_fg == 1){
		n_ref1 = 0.0;
		n_ref2 = 0.0;
		n_ref3 = -1.0*dz/fabs(dz);
		//printf("n_ref3 = %lg \n",n_ref3);
    }
    else{
		n_ref1 = -1.0 * dx;
		n_ref2 = -1.0 * dy;
		n_ref3 = -1.0 * dz;
		//printf("else\n");
	}

	ds1 = -1.0 * dx;
	ds2 = -1.0 * dy;
	ds3 = -1.0 * dz;

	theta = fabs(acos(ds1*n_ref1 + ds2*n_ref2 + ds3*n_ref3));

	if(theta>(PI/2) || theta<0){
		printf("dx = %g,dy = %lg, dz = %lg\n",dx,dy,dz);
		printf("du = %lg\n",sqrt(dx*dx + dy*dy + dz*dz));
		printf("theta=%lg\n",theta);
		printf("Len=%lg\n\n",Len);
	}

	duu=(ds1*n_ref1 + ds2*n_ref2 + ds3*n_ref3);

	dx = 2.0 * (ds1 * n_ref1 + ds2 * n_ref2 + ds3 * n_ref3) * n_ref1-ds1;
	dy = 2.0 * (ds1 * n_ref1 + ds2 * n_ref2 + ds3 * n_ref3) * n_ref2-ds2;
	dz = 2.0 * (ds1 * n_ref1 + ds2 * n_ref2 + ds3 * n_ref3) * n_ref3-ds3;

	refl = ref_out[(long)((fabs(duu * 1000)) + 0.5)];    /* calculate reflection */
	trans_ratio = 1.0 - refl;

	if (refl == 1)
	inr = 1;                              /* total internal reflection */
	else{
		temp = 1-(1-duu*duu) * ref_ind*ref_ind;
		if (temp < 0.0){                         /* can happen due to foundoff */
			inr = 1;
		}
		else{
			sweight = cweight * (1 - refl);
			cweight = cweight * refl;
			inr =0;
		}
	}
}


static void Output_no_ref(){
	/* calculate direction cosines relevant to tangent plane */
	sweight = cweight;
	cweight = 0.0;
	dx = -dx; //for acos in Store_data()
	dy = -dy;
	dz = -dz;
}


static void Source_pos(){
	double rand1,rand2,radius,theta;

	rand1 = GenRand();
	radius = source_minradius + ((source_maxradius - source_minradius) * rand1);
	rand2 = GenRand();
	theta = 2 * PI * rand2;
	if(source_dx == 1 || source_dx == -1){
		yold = yold + radius * cos(theta);
		zold = zold + radius * sin(theta);
	}
	else if(source_dy == 1 || source_dy == -1){
		xold = xold - radius * cos(theta);
		zold = zold + radius * sin(theta);
	}
	else if(source_dz == 1 || source_dz == -1){
		xold = xold + radius * cos(theta);
		yold = yold + radius * sin(theta);
	}
}


static void Source_dir(){
	/* calculate  direction vector at source (isotropic source)*/
	long ir1, ir2;
	double norm;

	ir1 = (long)(GenRand()*((double)TABLEN-1.0));
	ir2 = (long)(GenRand()*((double)TABLEN-1.0));

	if(fabs(source_dz) > COSZERO) {   /* normal launch */
		dx = source_sinthe[ir1]*source_cospsi[ir2];
		dy = source_sinthe[ir1]*source_sinpsi[ir2];
		if(source_dz >= 0){
			dz = source_costhe[ir1];
		}
		else {
			dz = - source_costhe[ir1];
		}
	}
	else {
		double temp = sqrt(1.0 - source_dz*source_dz);
		dx = source_sinthe[ir1] * (source_dx * source_dz * source_cospsi[ir2] - source_dy * source_sinpsi[ir2]) / temp + source_dx * source_costhe[ir1];
		dy = source_sinthe[ir1] * (source_dy * source_dz * source_cospsi[ir2] + source_dx * source_sinpsi[ir2]) / temp + source_dy * source_costhe[ir1];
		dz = -source_sinthe[ir1] * source_cospsi[ir2] * temp + source_dz * source_costhe[ir1];
	}
	norm = sqrt(dx*dx + dy*dy + dz*dz);
	dx /= norm;
	dy /= norm;
	dz /= norm;
}


/*[20180209�ǈ�]�Ǝˎ��Ԃ̈ʒu�̌���*/
static void Source_time(){
	double rand_time;

	rand_time = GenRand();
	if	(	0.00000000000000 	<=	rand_time	&&		rand_time	<=	1.00000000 	)	{	source_time	=	0	;	}




}


static void New_len(){
	long rand;

	rand = (long)(GenRand()*((double)TABLEN-2.0));
	Len = stepsize_table[rand] / mut[layer_new];
	if(Len == 0){
		printf("rand = %ld\nLen = %f\n",rand,stepsize_table[rand]);
		Len = V0;
	}
}


static void New_dir(){
	/* calculate new direction vector (isotropic scattering)*/
	double dr1, dr2, dr3;     /* relative new direction cosines */
	double dv1, dt3;
	long ir1, ir2;

	ir1 = (long)(GenRand()*((double)TABLEN-1.0));
	ir2 = (long)(GenRand()*((double)TABLEN-1.0));

	if (!strcmp(angfile,"isotropic")){
		dx = sinthe[ir1]*cospsi[ir2];
		dy = sinthe[ir1]*sinpsi[ir2];
		dz = costhe[ir1];
	}
	else{
		dr1 = sinthe[ir1]*cospsi[ir2];
		dr2 = sinthe[ir1]*sinpsi[ir2];
		dr3 = costhe[ir1];
		dt3 = 1-dz*dz;
		if (dt3>0){
			dt3 = sqrt(dt3);
			dv1 = (dx*dz*dr1-dy*dr2)/dt3+dx*dr3;
			dy = (dy*dz*dr1+dx*dr2)/dt3+dy*dr3;
			dz = dz*dr3-dt3*dr1;
			dx = dv1;
			return;
		}
		if (dz>0){            /* same direction, theta=0 */
			dx = dr1;
			dy = dr2;
			dz = dr3;
			return;
		}
		dx = -dr1;             /* oposit direction theta=180 deg.*/
		dy = -dr2;
		dz = -dr3;
	}
}


static void Check_layer(double z){
	int layer_num,layer_min;
	layer_min = 0;
	if(z < 0){
		printf("ERROR:z = %f\n", z);
		exit(1);
	}
	for(layer_num = 0; layer_num < MAX_LAYER; layer_num++){
		if(z >= layer_min && z < (layer_min + thickness[layer_num])){
			layer_new = layer_num;
			break;
		}
		else if(layer_num == MAX_LAYER - 1){
			layer_new = layer_num;
		}
		else{
			layer_min += thickness[layer_num];
		}
	}
}


static void Restore_voxelpath(){
	int n,restore_fg;
	double xoldtemp, yoldtemp, zoldtemp;
	double L,L_temp,Lx,Ly,Lz;

	restore_fg = 0;
	L_temp = 0.0;
	xoldtemp = xold;
	yoldtemp = yold;
	zoldtemp = zold;
	if(xoldtemp < 0){	vox = (int)(xoldtemp) - 1;}
	else{				vox = (int)(xoldtemp);}
	if(yoldtemp < 0){	voy = (int)(yoldtemp) - 1;}
	else{				voy = (int)(yoldtemp);}
	voz = (int)(zoldtemp);
	n = 0;
	L = Len;

	do{
        printf("hello\n");

		Check_layer(zoldtemp);
		layer_old = layer_new;

		refx_fg = 0;
		refy_fg = 0;
		refz_fg = 0;

		if(dx > 0){
			Lx = ((double)vox+1.0-xoldtemp)/dx;
		}
		else if (dx < 0) {
			Lx = ((double)vox-xoldtemp)/dx;
		}
		else if(dx==0) {
			Lx=1.0/0.0;
		}

		if(dy > 0){
			Ly= ((double)voy+1.0-yoldtemp)/dy;
		}
	    else if(dy < 0){
			Ly= ((double)voy-yoldtemp)/dy;
		}
		else if(dy==0){
			Ly=1.0/0.0;
		}

		if(dz > 0){
			Lz= ((double)voz+1.0-zoldtemp)/dz;
		}
		else if (dz < 0) {
			Lz= ((double)voz-zoldtemp)/dz;
		}
		else if(dz == 0) {
			Lz=1.0/0.0;
		}

		if ((Lx < Ly) && (Lx < Lz) && (Lx <= L)){
			L_temp += Lx;
			path[layer_old] += Lx;
			if((PD_MIN_LINE <= voy) && (voy <= PD_MAX_LINE) && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){
				path_for_ssp[voz][vox+MAX_X][voy - PD_MIN_LINE] += Lx;

			if(voy == (int)source_y && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){

				path_y[layer_old] += Lx;
			}
				}


			L = L - Lx;
			refx_fg = 1;
			xoldtemp = xoldtemp+Lx*dx;
			yoldtemp = yoldtemp+Lx*dy;
			zoldtemp = zoldtemp+Lx*dz;
			if (dx>0){
			  vox +=1;
			}
			else if(dx<0){
				vox -=1;
			}
		}
		else if((Ly < Lx) && (Ly < Lz) && (Ly <= L)){
			L_temp += Ly;
			path[layer_old] += Ly;
			if((PD_MIN_LINE <= voy) && (voy <= PD_MAX_LINE) && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){
				path_for_ssp[voz][vox+MAX_X][voy - PD_MIN_LINE] += Ly;

			if(voy == (int)source_y && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){
				path_y[layer_old] += Ly;
			}
			}
			L = L - Ly;
			refy_fg = 1;
			xoldtemp = xoldtemp+Ly*dx;
			yoldtemp = yoldtemp+Ly*dy;
			zoldtemp = zoldtemp+Ly*dz;
			if (dy > 0){
				voy +=1;
			}
			else if (dy < 0){
				voy -=1;
			}
		}
		else if((Lz < Lx) && (Lz < Ly) && (Lz <= L)){
			L_temp += Lz;
			path[layer_old] += Lz;
			if((PD_MIN_LINE <= voy) && (voy <= PD_MAX_LINE) && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){
				path_for_ssp[voz][vox+MAX_X][voy - PD_MIN_LINE] += Lz;

			if(voy == (int)source_y && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){
				path_y[layer_old] += Lz;
			}
			}
			L = L - Lz;
			refz_fg = 1;
			xoldtemp = xoldtemp+Lz*dx;
			yoldtemp = yoldtemp+Lz*dy;
			zoldtemp = zoldtemp+Lz*dz;
			if (dz>0){
			    voz +=1;
			}
			else if(dz < 0){
			    voz -=1;
			}
		}
		else if((Lx == Ly) && (Lx < Lz) && (Lx <= L)){
			L_temp += Lx;
			path[layer_old] += Lx;
			if((PD_MIN_LINE <= voy) && (voy <= PD_MAX_LINE) && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){
				path_for_ssp[voz][vox+MAX_X][voy - PD_MIN_LINE] += Lx;

			if(voy == (int)source_y && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){

				path_y[layer_old] += Lx;
			}
			}
			L = L - Lx;
			refx_fg = 1;
			refy_fg = 1;
			xoldtemp = xoldtemp+Lx*dx;
			yoldtemp = yoldtemp+Lx*dy;
			zoldtemp = zoldtemp+Lx*dz;
			if (dx>0){
				vox +=1;
			}
			else if(dx<0){
				vox -=1;
			}
			if (dy>0){
				voy +=1;
			}
			else if(dy<0){
				voy -=1;
			}
		}
		else if((Lx == Lz) && (Lx < Ly) && (Lx <= L)){
			L_temp += Lx;
			path[layer_old] += Lx;
			if((PD_MIN_LINE <= voy) && (voy <= PD_MAX_LINE) && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){
				path_for_ssp[voz][vox+MAX_X][voy - PD_MIN_LINE] += Lx;

			if(voy == (int)source_y && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){
				path_y[layer_old] += Lx;
			}
			}
			L = L - Lx;
			refx_fg = 1;
			refz_fg = 1;
			xoldtemp = xoldtemp+Lx*dx;
			yoldtemp = yoldtemp+Lx*dy;
			zoldtemp = zoldtemp+Lx*dz;
			if (dx>0){
				vox +=1;
			}
			else if(dx<0){
				vox -=1;
			}
			if (dz>0){
				voz +=1;
			}
			else if(dz<0){
				voz -=1;
			}
		}
		else if((Ly == Lz) && (Ly < Lx) && (Ly <= L)){
			L_temp += Ly;
			path[layer_old] += Ly;
			if((PD_MIN_LINE <= voy) && (voy <= PD_MAX_LINE) && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){
				path_for_ssp[voz][vox+MAX_X][voy - PD_MIN_LINE] += Ly;

			if(voy == (int)source_y && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){

				path_y[layer_old] += Ly;
			}
			}
			L = L - Ly;
			refy_fg = 1;
			refz_fg = 1;
			xoldtemp = xoldtemp+Ly*dx;
			yoldtemp = yoldtemp+Ly*dy;
			zoldtemp = zoldtemp+Ly*dz;
			if (dy>0){
				voy +=1;
			}
			else if(dy<0){
				voy -=1;
			}
			if (dz>0){
				voz +=1;
			}
			else if(dz<0){
				voz -=1;
			}
		}
		else if((Lx == Ly) && (Ly == Lz) && (Lz == Lx) && (Lx <= L)){
			L_temp += Lx;
			path[layer_old] += Lx;
			if((PD_MIN_LINE <= voy) && (voy <= PD_MAX_LINE) && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){
				path_for_ssp[voz][vox+MAX_X][voy - PD_MIN_LINE] += Lx;

			if(voy == (int)source_y && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){

				path_y[layer_old] += Lx;
			}
			}
			L = L - Lx;
			refx_fg = 1;
			refy_fg = 1;
			refz_fg = 1;
			xoldtemp = xoldtemp+Lx*dx;
			yoldtemp = yoldtemp+Lx*dy;
			zoldtemp = zoldtemp+Lx*dz;

			if (dx>0){
				vox +=1;
			}
			else if(dx<0){
				vox -=1;
			}
			if (dy>0){
				voy +=1;
			}
			else if(dy < 0){
				voy -=1;
			}
			if (dz>0){
				voz +=1;
			}
			else if(dz < 0){
				voz -=1;
			}
		}
		else if((L < Lx) && (L < Ly) && (L < Lz)){
			L_temp += L;
			path[layer_old] += L;
			if((PD_MIN_LINE <= voy) && (voy <= PD_MAX_LINE) && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){
				path_for_ssp[voz][vox+MAX_X][voy - PD_MIN_LINE] += L;

			if(voy == (int)source_y && voz >= 0 && voz < MAX_Z && vox >= -MAX_X && vox < MAX_X+1){

				path_y[layer_old] += L;
			}
			}
			L=0.0;
			restore_fg = 1;
		}
		else {
			printf("ERROR\n");
		}

		if(voz < 0){
			L_temp +=L;
			restore_fg = 1;
			zoldtemp = 0;
			cross_x = xoldtemp;
			cross_y = yoldtemp;
			cross_z = zoldtemp;
		}
		if(zoldtemp < 0){
			printf("Restore_ERROR\nzoldtemp = %f,voz = %d\n",zoldtemp,voz);
			exit(1);
		}
		Check_layer(zoldtemp);
		if(layer_new != layer_old){
			L = L*mut[layer_old]/mut[layer_new];
		}
	}while(restore_fg == 0);

	Len = L_temp;
}


static void Restore_data(){
	/* correct total pathlength and find rigth X- and Y coordinates */
	Lsub = sqrt((xnew - cross_x)*(xnew -cross_x) + (ynew - cross_y)*(ynew -cross_y) + (znew - cross_z)*(znew - cross_z));
	allpath -= Lsub;
	Len = Lsub;

	if (Len == 0){
		Len = (double)V0;
		Len = -(log(Len)/mut[0])/SCALE;
	}

	xnew= cross_x;
	ynew= cross_y;
	znew= cross_z;
}

static void Store_data(){
	int		i,j,x,y,z;
	long	time,time_SSP;
	double	d1,d2,inner_product,theta;

	time = labs((long)(((allpath*ref_ind)*SCALE / VC + source_time )/ DT));
	time_SSP = labs((long)(((allpath*ref_ind)*SCALE / VC + source_time) / DT_SSP));

/*	if (source_time == 1){
		printf("sorce_time_ok\n");

	}
*/
	for(i = 0; i < MAX_DET; i++){
		inner_product = detector_dx[i] * dx + detector_dy[i] * dy + detector_dz[i] * dz;
		theta = acos(inner_product);

		if(detector_dz[i] == 1 || detector_dz[i] == -1){
			d1 = (detector_x[i] - xnew) * (detector_x[i] - xnew);
			d2 = (detector_y[i] - ynew) * (detector_y[i] - ynew);
		}

		if(sqrt(d1+d2) >= (double)detector_minradius/SCALE && sqrt(d1+d2) <= (double)detector_maxradius/SCALE && theta <= asin(detector_NA)){
			phot_detected++;
			for (j=0; j<MAX_LAYER; j++){
				ppath[i][j] += sweight*path[j]*SCALE;

				if((time >= 0) && (time < MT)){
					tppath[i][time][j] += sweight*path[j]*SCALE;
					tppath_y[i][time][j] += sweight*path_y[j] * SCALE;
				}
			}
			for(z=0;z<MAX_Z;z++){
				for(x=0;x<MAX_X*2;x++){
									for(y=0;y<=MAX_Y*2;y++){
					SSP[i][z][x][y] += sweight * path_for_ssp[z][x][y] * SCALE;
									}
				}
			}
			if((time_SSP >= 0) && (time_SSP < MT_SSP)){
				intensity_for_ssp[i][time_SSP] += sweight;
				for(z=0;z<MAX_Z;z++){
					for(x=0;x<MAX_X*2;x++){
						for(y=0;y<=MAX_Y*2;y++){
						TSSP[i][time_SSP][z][x][y] += sweight * path_for_ssp[z][x][y] * SCALE;
						}
					}
				}
			}

			if((time >= 0) && (time <MT)){
				intensity[i][time] += sweight;
			}
			return;
		}
	}
}


static void Fnstop(){
	int i;
	unsigned long long int sec;
	int month,day,hour,min;

	time_end = time(NULL);
	phot_end = phot_in;
	sec = (double)(time_end - time_start)/(double)(phot_end - phot_start)*1000000000.0;
	if((phot_end - phot_start) == 0){
		sec = 0;
	}
	month = sec / (60*60*24*30);
	sec = sec % (60*60*24*30);
	day = sec / (60*60*24);
	sec = sec % (60*60*24);
	hour = sec / (60*60);
	sec = sec % (60*60);
	min = sec / 60;
	sec = sec % 60;

	fp_com = fopen(comfile, "r");
	if (fp_com == NULL){
		fprintf(stderr,"comfile can not open\n");
		exit(1);
	}

	fscanf(fp_com, "%c%*[^\n]", &ch);
	getc(fp_com);
	fclose(fp_com);

	fp_com = fopen(comfile, "w");
	if (fp_com == NULL){
		fprintf(stderr,"comfile can not open\n");
		exit(1);
	}

	fprintf(fp_com, "<This file is created by 'mc_fnirs_slab.c'>\n");
	fprintf(fp_com, "number of input photons \t %12llu\n", phot_in);
	fprintf(fp_com, "number of output photons \t %12llu\n", phot_out);
	fprintf(fp_com, "number of detect photons \t %12llu\n", phot_detected);
	fprintf(fp_com, "number of over time photon \t %12llu\n", phot_overtime);
	fprintf(fp_com, "number of over path photon \t %12llu\n", phot_overpath);
	fprintf(fp_com, "number of over scatter photon \t %12llu\n", phot_overscat);
	fprintf(fp_com, "estimate time of 1B photon = %dM, %dD, %dH, %dM, %lluS\n\n",month,day,hour,min,sec);

	fprintf(fp_com, "Layer Num\tThickness\tmus\tmua\n");
	for (i=0;i<MAX_LAYER;i++){
		fprintf(fp_com, "%d\t\t%d\t\t%6.4f\t%6.4f\n", i, thickness[i], mus[i], mua[i]);
	}

	fprintf(fp_com, "\nPhase function\t: %s (g = %.1lf)\n", angfile,g);
	fprintf(fp_com, "Source NA\t: %.1f\n", source_NA);
	fprintf(fp_com, "Source radius\t: Min = %.1f\tMax = %.1f\n", source_minradius, source_maxradius);
	fprintf(fp_com, "Source\t\t: (x = %.1f, y = %.1f, z = %.1f)\n", source_x,source_y,source_z);
	fprintf(fp_com, "Detector NA\t: %.1f\n", detector_NA);
	fprintf(fp_com, "Detector radius\t: Min = %.1f\tMax = %.1f\n", detector_minradius,detector_maxradius);
	for(i = 0; i < MAX_DET; i++){
		fprintf(fp_com, "Detector[%d]\t: (x = %.1f, y = %.1f, z = %.1f)\n", i, detector_x[i], detector_y[i], detector_z[i]);
	}

	if ((ch=='y')||(ch=='Y'))
	{
		stop_fg = TRUE;
		fprintf(fp_com, "stop request");
	}
	else{
		stop_fg = FALSE;
	}

	fclose(fp_com);
}
