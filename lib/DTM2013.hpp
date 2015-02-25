/*
 * Header file for calling DTM2012 from C++
 */

#ifndef DTM2013_H_
#define DTM2013_H_



//Custom type required to use "dtm_wrapper"
struct dtm_date {
   int type_flag;
   double mjd2000;
   int day;
   int month;
   int year;
   int hour;
   int minute;
   double second;
};

// DTM2012 subroutines
extern "C" {
void dtm_wrapper (struct dtm_date *in_date, float *alti, float *alat, float *xlon, float *tz, float *tinf, float *tp120, float *ro, float *ro_unc, float *d, float *wmm);
void load_config (char filePath[200]);
void dtm2013 (float *day,float *f, float *fbar, float *akp, float *alti, float *hl, float *alat, float *xlon,float *tz, float *tinf, float *tp120, float *ro, float *d, float* wmm);
void density_uncertainty (float *alt, float *lat, float *lst, float *flux, float *kp, float *unc);
void P_ReadDTM13(void);
}

#endif /* DTM2013_H_ */
