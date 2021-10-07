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
#define STRSIZE		50
typedef char string[STRSIZE];

static FILE		*fp_base, *fp_model, *fp_note, *fp_data, *fp_pd, *fp_ssp;
static string	baseconf = "settings.conf";
static string	modelconf = "model.conf";
static string	note = "temporary_note.txt";
static string	datafile = "binary.data";
static string	pdfile = "binary.pd";
static string	sspfile = "binary.ssp";
/* ========================================================================= */

/* global variables to be altered ========================================== */
#define SCALE		1.0		/* size of pixel*/

#define MAX_X		60		/* data save x (-MAX_X < MAX_X)*/
#define MAX_Y		2		/* data save x (-MAX_X < MAX_X)*/
#define MAX_Z		50		/* data save z (0 < MAX_Z)*/

#define PD_MAX_LINE	5
#define PD_MIN_LINE	1		/*y�������̋��߂�PD�̍��W�͈́i-80�`80�j�i�͈͂�10�ȓ��j*/

#define MAX_DET		1
#define MAX_LAYER	5

#define MT			400		/* number of division for TPSF */
#define DT			10
#define MT_PD		80		/* number of division for TR photon density*/
#define DT_PD		50
#define MT_SSP		16		/* number of division for TR spatial */
							/*sencitivity profile*/
#define DT_SSP		250
#define INT_CUT		0.0001	/* cutoff intensity for cweight */
#define MAX_SCATTER	100000
#define TABLEN		10000	/* length of lookup table for new direction cos */
/* ========================================================================= */

/* global variables as constants =========================================== */
#define VC			0.3				/* speed of light in vacuum (mm/psec) */
#define REF_IND		1.4

#define PI			3.14159265358979324
#define DEG_RAD		(PI / 180.0)
#define RAD_DEG		(180.0 / PI)
#define COSZERO		(1.0 - 1.0E-12)			/* cosine of about 1e-6 rad */

#define V0			1.77635683940025e-15	/* =2**-49 (min of randnumber) */

#define TRUE		1
#define FALSE		0

// extern time_t time(time_t *);
#define N				25
#define M				7
static unsigned long	x[N] = {
	0x95f24dab, 0x0b685215, 0xe76ccae7, 0xaf3ec239, 0x715fad23,
	0x24a590ad, 0x69e4b5ef, 0xbf456141, 0x96bc1b7b, 0xa7bdf825,
	0xc1de75b7, 0x8858a9c9, 0x2da87693, 0xb657f9dd, 0xffdc8a9f,
	0x8121da71, 0x8b823ecb, 0x885d05f5, 0x4e20cd47, 0x5a9ad5d9,
	0x512c0c03, 0xea857ccd, 0x4cc1d30f, 0x8891a8a1, 0xa6b7aadb
};
static unsigned long	mag[2] = {0x0, 0x8ebfd028};
static int k = 0;
/* ========================================================================= */

/* variables to read from file ============================================= */
static char						is_newfile;	/* if new file (y/n) */
static unsigned long long int	phot_input;	/* total number of input photons */
static double					g;			/* mean cosine of phase function */
static double					maxpath; 	/* max pathlength for cutoff */

static double	source_NA;
static double	source_minradius;
static double	source_maxradius;
static double	source_x;
static double	source_y;
static double	source_z;
static double	source_dx;
static double	source_dy;
static double	source_dz;

static double	detector_NA;
static double	detector_minradius;
static double	detector_maxradius;
static double	detector_x[MAX_DET];
static double	detector_y[MAX_DET];
static double	detector_z[MAX_DET];
static double	detector_dx[MAX_DET];
static double	detector_dy[MAX_DET];
static double	detector_dz[MAX_DET];

static int		thickness[MAX_LAYER];
static double	mus[MAX_LAYER];
static double	mua[MAX_LAYER];
static double	mut[MAX_LAYER];
/* ========================================================================= */

/* variables for recording ================================================= */
static unsigned long long int	phot_start, phot_end;
static time_t					time_start, time_end;

static unsigned long long int	phot_in;		/* input photons */
static unsigned long long int	phot_out;		/* output photons */
static unsigned long long int	phot_detected;	/* detected photons */
static unsigned long long int	phot_overtime;	/* over-time photons */
static unsigned long long int	phot_overpath;	/* over-path photons */
static unsigned long long int	phot_overscat;	/* over-scattering photons */

static double	PD[MAX_Z][MAX_X*2+1];
static double	TPD[MT_PD][MAX_Z][MAX_X * 2 + 1];
static double	STPD[MT_PD];
static double	SSP[MAX_DET][MAX_Z][MAX_X*2+1][MAX_Y*2+1];
static double	TSSP[MAX_DET][MT_SSP][MAX_Z][MAX_X*2+1][MAX_Y*2+1];

static double	path[MAX_LAYER];
static double	path_y[MAX_LAYER];
static double	path_for_ssp[MAX_Z][MAX_X*2+1][MAX_Y * 2 + 1];
static double	ppath[MAX_DET][MAX_LAYER];
static double	tppath[MAX_DET][MT][MAX_LAYER];
static double	tppath_y[MAX_DET][MT][MAX_LAYER];
static double	intensity[MAX_DET][MT];
static double	intensity_for_ssp[MAX_DET][MT_SSP];
/* ========================================================================= */

/* subfunction declaration ================================================= */
static void LoadSettings(); /* loads values from settings file */
static void InitData(); /* initializes variables if file new */
static void LoadData(); /* loads variables if result partially recorded */
static void SaveData(); /* saves variables to binary data */
static void SetSeed();
static double GenRand();
/* ========================================================================= */

/* dummy variables ========================================================= */
long	i, j;
/* ========================================================================= */


int main(void){
	LoadSettings();

	if(is_newfile == 'n'){
		LoadData();
	}
	else{
		InitData();
	}

	// New_dir_table();
	// Source_dir_table();
	// New_len_table();
	// Ref_table();
	SetSeed();
	printf("%lf\n", GenRand());
	printf("%lf\n", GenRand());

	time_start = time(NULL);
	// Main_mc();
	// Fnstop();
	// SaveData();

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
