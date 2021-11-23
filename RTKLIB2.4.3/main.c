/*------------------------------------------------------------------------------
* rnx2rtkp.c : read rinex obs/nav files and compute receiver positions
*
*          Copyright (C) 2007-2009 by T.TAKASU, All rights reserved.
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:55:16 $
* history : 2007/01/16  1.0 new
*           2007/03/15  1.1 add library mode
*           2007/05/08  1.2 separate from postpos.c
*           2009/01/20  1.3 support rtklib 2.2.0 api
*           2009/12/12  1.4 support glonass
*                           add option -h, -a, -l, -x
*           2010/01/28  1.5 add option -k
*           2010/08/12  1.6 add option -y implementation (2.4.0_p1)
*           2014/01/27  1.7 fix bug on default output time format
*-----------------------------------------------------------------------------*/
#include <stdarg.h>
#include "rtklib.h"

static const char rcsid[] = "$Id: rnx2rtkp.c,v 1.1 2008/07/17 21:55:16 ttaka Exp $";

#define PROGNAME    "rnx2rtkp"          /* program name */
#define MAXFILE     8                   /* max number of input files */

/* help text -----------------------------------------------------------------*/
static const char *help[] = {
"",
" usage: rnx2rtkp [option]... file file [...]",
"",
" Read RINEX OBS/NAV/GNAV/HNAV/CLK, SP3, SBAS message log files and ccompute ",
" receiver (rover) positions and output position solutions.",
" The first RINEX OBS file shall contain receiver (rover) observations. For the",
" relative mode, the second RINEX OBS file shall contain reference",
" (base station) receiver observations. At least one RINEX NAV/GNAV/HNAV",
" file shall be included in input files. To use SP3 precise ephemeris, specify",
" the path in the files. The extension of the SP3 file shall be .sp3 or .eph.",
" All of the input file paths can include wild-cards (*). To avoid command",
" line deployment of wild-cards, use \"...\" for paths with wild-cards.",
" Command line options are as follows ([]:default). With -k option, the",
" processing options are input from the configuration file. In this case,",
" command line options precede options in the configuration file.",
"",
" -?        print help",
" -k file   input options from configuration file [off]",
" -o file   set output file [stdout]",
" -ts ds ts start day/time (ds=y/m/d ts=h:m:s) [obs start time]",
" -te de te end day/time   (de=y/m/d te=h:m:s) [obs end time]",
" -ti tint  time interval (sec) [all]",
" -p mode   mode (0:single,1:dgps,2:kinematic,3:static,4:moving-base,",
"                 5:fixed,6:ppp-kinematic,7:ppp-static) [2]",
" -m mask   elevation mask angle (deg) [15]",
" -f freq   number of frequencies for relative mode (1:L1,2:L1+L2,3:L1+L2+L5) [2]",
" -v thres  validation threshold for integer ambiguity (0.0:no AR) [3.0]",
" -b        backward solutions [off]",
" -c        forward/backward combined solutions [off]",
" -i        instantaneous integer ambiguity resolution [off]",
" -h        fix and hold for integer ambiguity resolution [off]",
" -e        output x/y/z-ecef position [latitude/longitude/height]",
" -a        output e/n/u-baseline [latitude/longitude/height]",
" -n        output NMEA-0183 GGA sentence [off]",
" -g        output latitude/longitude in the form of ddd mm ss.ss' [ddd.ddd]",
" -t        output time in the form of yyyy/mm/dd hh:mm:ss.ss [sssss.ss]",
" -u        output time in utc [gpst]",
" -d col    number of decimals in time [3]",
" -s sep    field separator [' ']",
" -r x y z  reference (base) receiver ecef pos (m) [average of single pos]",
" -l lat lon hgt reference (base) receiver latitude/longitude/height (deg/m)",
" -y level  output soltion status (0:off,1:states,2:residuals) [0]",
" -x level  debug trace level (0:off) [0]"
};


/* print help ----------------------------------------------------------------*/
static void printhelp(void)
{
	int i;
	for (i = 0; i < (int)(sizeof(help) / sizeof(*help)); i++) fprintf(stderr, "%s\n", help[i]);
	exit(0);
}

/*
m:band
n:obs number
*/
void CalcStats(const sol_t *sol, int n, int m, double ref, double *ave, double *std, double *rms, double *cep) {
	double sum = 0, sumsq = 0;
	int i;
	if (n <= 0) {
		return;
	}

	for (int i = 0; i < n; i++) {
		sum += sol[i].rr[m];
		sumsq += sol[i].rr[m] * sol[i].rr[m];
	}
	*ave = sum / n;
	*std = n > 1 ? sqrt((sumsq - 2 * sum*(*ave) + *ave * *ave*n) / (n - 1)) : 0;
	*rms = sqrt((sumsq - 2 * sum*ref + ref * ref*n) / n);
	*cep = *rms / 1.2;
}
void calRMS(solbuf_t *solbuf, int n) {
	/*
	n=0 开始点为参考 1：结束点为参考 2：均值点为参考
	*/

	sol_t *sol = NULL;
	int m = 0;
	double ave[3] = { 0 }, std[3] = { 0 }, rms[3] = { 0 }, cep[3] = { 0 };
	double opos[3] = { 0 };
	if (n == 0)
	{
		if (!(sol = getsol(solbuf, 0) && sol->type != 0)) return;
		for (int i = 0; i < 3; i++) { opos[i] = sol->rr[i]; }

	}
	else if (n == 1) {
		if (!(sol = getsol(solbuf, solbuf[0].n - 1) && sol->type != 0)) { return; }
		for (int i = 0; i < 3; i++) { opos[i] = sol->rr[i]; }
	}
	else if (n == 2) {
		for (int i = 0; i < solbuf[0].n - 1; i++) {
			if (solbuf->data[0].type != 0) {
				continue;
			}for (int j = 0; j < 3; j++) {
				opos[j] += solbuf->data[i].rr[j];
			}
			m++;
		}
		if (m > 0)for (int i = 0; i < 3; i++) { opos[i] /= m; }
	}
	CalcStats(solbuf->data, solbuf[0].n, 0, opos[0], &ave[0], &std[0], &rms[0], &cep[0]);
	CalcStats(solbuf->data, solbuf[0].n, 1, opos[1], &ave[1], &std[1], &rms[1], &cep[1]);
	CalcStats(solbuf->data, solbuf[0].n, 2, opos[2], &ave[2], &std[2], &rms[2], &cep[2]);
}
/* rnx2rtkp main -------------------------------------------------------------*/
int main()
{
	prcopt_t prcopt = prcopt_default;
	solopt_t solopt = solopt_default;
	filopt_t filopt = { "" };
	gtime_t ts = { 0 }, te = { 0 };
	//double st[] = { 2021, 4, 28, 5, 51, 22 }, se[] = { 2021, 4, 28, 6, 46, 48 };

	//char* infile[] = { {"D:\\作业\\手机RTK\\实验数据\\华为mate30\\4月28日 大东门\\Geo++118F..21o"},{"D:\\作业\\手机RTK\\实验数据\\RTK\\4月28日 大东门\\DDM_0428_1.21p"},{"D:\\作业\\手机RTK\\实验数据\\RTK/4月28日 大东门\\gfz21553.sp3"} };
	//char* outfile = { "D:\\作业\\手机RTK\\实验数据\\华为mate30\\4月28日 大东门\\Geo++118F2..pos" };
	char* infile[] = { {"D:\\作业\\手机RTK\\实验数据\\RTK\\4月28日 大东门\\DDM_0428_1.21o"},{"D:\\作业\\手机RTK\\实验数据\\RTK\\4月28日 大东门\\DDM_0428_1.21p"},{"D:\\作业\\手机RTK\\实验数据\\RTK\\4月28日 大东门\\gfz21553.sp3"} };
	char* outfile = { "D:\\作业\\手机RTK\\实验数据\\RTK\\4月28日 大东门\\DDM_0428_1.pos" };
	double tint = 0.0, pos[3];
	int i, j, n = 3, ret;
	//te = epoch2time(se);
	//ts = epoch2time(st);
	prcopt.mode = PMODE_PPP_STATIC;
	prcopt.navsys = SYS_ALL;
	prcopt.ionoopt = IONOOPT_BRDC;
	prcopt.tropopt = TROPOPT_SAAS;
	prcopt.sateph = EPHOPT_BRDC;

	solopt.posf = SOLF_LLH;
	strcpy(solopt.prog, "iRefer");
	solopt.sstat = 2;
	solopt.trace = 5;
	ret = postpos(ts, te, tint, 0.0, &prcopt, &solopt, &filopt, infile, n, outfile, "", "");

	solbuf_t solbuf = { 0 };
	sol_t *sol;
	char* files[] = { "D:/作业/手机RTK/实验数据/华为mate30/4月28日 大东门/Geo++118F2..pos" };

	ret = readsolt(files, 1, ts, te, tint, 0, &solbuf);
	calRMS(&solbuf, 2);

	if (!ret) fprintf(stderr, "%40s\r", "");
	return ret;
}