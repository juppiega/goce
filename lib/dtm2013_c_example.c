#include<stdio.h>
#include<string.h>
#include"DTM2013.h"

int main() {

	//Path to configuration file
	char configFile[200];
	strncpy(configFile,"config.cfg",200);

	//DTM wrapper inputs
    struct dtm_date today; //struct provided by DTM2012.h
    float alti;
    float alat;
    float xlon;

    //DTM wrapper outputs
    float ro;
    float tinf;
    float tz;
    float tp120;
    float wmm;
    float ro_unc;
    float d[6];

    //Fill in the universal time in calendar date format
    today.type_flag=2;
    today.day=3;
    today.month=9;
    today.year=2011;
    today.hour=8;
    today.minute=32;
    today.second=34.23;

    //Fill in the other inputs
    alti=200.0;
    alat=40.0;
    xlon=-4.0;

    //load configuration
    load_config(configFile);

    //Call dtm wrapper
    dtm_wrapper (&today, &alti, &alat, &xlon, &tz, &tinf, &tp120, &ro, &ro_unc, d, &wmm);

    printf("\n");
    printf("-----------CALL TO DTM2012-----------\n");
    printf("--------------C WRAPPER--------------\n");
    printf("\n");


    printf("inputs\n");
    printf("   Date:             %f (%i/%i/%i) at %i : %i : %f\n", today.mjd2000, today.day, today.month,today.year, today.hour, today.minute, today.second);
    printf("\n");
    printf("\n");
    printf("   latitude:            %f\n", alat);
    printf("   longitude:           %f\n", xlon);
    printf("   altitude:            %f\n", alti);
    printf("\n");
    printf("\n");


    printf("\n");
    printf("outputs\n");
    printf("   Temp at altitude : %f\n", tz);
    printf("   exospheric tmp :   %f\n", tinf);
    printf("   atomic hydrogen :  %g\n", d[0]);
    printf("   helium :           %g\n", d[1]);
    printf("   atomic oxygen :    %g\n", d[2]);
    printf("   molecular nitro :  %g\n", d[3]);
    printf("   molecular oxygen : %g\n", d[4]);
    printf("   atomic nitrogen :  %g\n", d[5]);
    printf("   density (g/cm^3) : %g   (+-  %f %% )\n", ro, ro_unc);
    printf("   mean molec mass :  %g\n", wmm);


	return 0;
}
