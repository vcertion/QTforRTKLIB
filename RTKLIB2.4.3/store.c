#include "rtklib.h"
double P[MAXSAT];
float D[MAXSAT];
unsigned char sat[MAXSAT];
time_t timet[MAXSAT];
nav_t nav_s[MAXSAT];
int number = 0, l_s, haveD = 0;
FILE *fp, *fp2;
obsd_t *obs_s;
double *rs_s, *vare_s, *svh_s, *dts_s;
char posMode[1], inoCore[1], troCore[1], sateCore[1];
char* infile[4];
extern void set(obsd_t *obs, nav_t *nav, int n, int l) {
	for (int i = 0; i < n; i++) {
		sat[i] = obs[i].sat;
		timet[i] = obs[i].time.time;
		P[i] = obs[i].P[l];
		D[i] = obs[i].D[l];
		l_s = l;
	}
	number = n;
}
extern int getNumber_s() {
	return number;
}
extern double getP(int n) {
	return P[n];
}
extern unsigned char getSat(int n) {
	return sat[n];
}
extern time_t getTime_s(int n) {
	return timet[n];
}
extern float getD(int n) {
	return D[n];
}
extern int getl() {
	return l_s;
}
extern void setHaveD(int n) {
	haveD = n;
}
extern int getHaveD() {
	return haveD;
}
extern void setobsOutFile(FILE *file) {
	fp = file;
}

extern FILE* getobsOutFile() {
	return fp;
}

extern void setRaimData(obsd_t *obs, int n, double *rs, double *dts, double *vare, double *svh, int j, int k) {
	obs_s = (obsd_t *)malloc(sizeof(obsd_t)*n);
	//obs_s[k] = obs;
	matcpy(rs + 6 * k, rs + 6 * j, 6, 1);
	matcpy(dts + 2 * k, dts + 2 * j, 2, 1);
	//vare_s[k] = vare;
}
extern void setPositingMode(char pm){
	posMode[0] = pm;
}
extern char getPositiongMode() {
	return posMode[0];
}
extern void setInoCore(char ic) {
	inoCore[0] = ic;
}
extern char getInoCore() {
	return inoCore[0];
}
extern void setTroCore(char tc) {
	troCore[0] = tc;
}
extern char getTroCore() {
	return troCore[0];
}
extern void setSateCore(char sc) {
	sateCore[0] = sc;
}
extern char getSateCore() {
	return sateCore[0];
}
extern void setInFile(char* in[4]) {
	for (int i = 0; i < 4; i++) {
		infile[i] = in[i];
	}
}
extern char* getInFile() {
	return infile;
}
