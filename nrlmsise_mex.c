#include "mex.h"
#include "nrlmsise-00.h"
#include <math.h>

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{
    // Declare the variables first because this is c
    double time, doy_in, year_in, sec_in, alt_in, lat_in, lon_in;
    double f107A_in, f107_in, ap_in, lst_in, ApDaily_in;
    struct nrlmsise_output output;
    struct nrlmsise_input input;
    struct ap_array apVals;
    struct nrlmsise_flags flags;
    int i;


    if(!(nrhs == 15 || nrhs == 9))
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","15 or 9 inputs required.");
    }

    // Get input

    doy_in           = mxGetScalar(prhs[0]);
    sec_in           = mxGetScalar(prhs[1]);
    alt_in           = mxGetScalar(prhs[2]);
    lat_in           = mxGetScalar(prhs[3]);
    lon_in           = mxGetScalar(prhs[4]);
    lst_in           = mxGetScalar(prhs[5]);
    f107A_in         = mxGetScalar(prhs[6]);
    f107_in          = mxGetScalar(prhs[7]);
    ApDaily_in       = mxGetScalar(prhs[8]);
      
    i;
    flags.switches[0] = 0;
    for(i = 1; i < 24; i++)
    {
        flags.switches[i] = 1;
    }
    
    if(nrhs == 15)
    {
        flags.switches[9] = -1;
        apVals.a[0]       = mxGetScalar(prhs[8]);
        apVals.a[1]       = mxGetScalar(prhs[9]);
        apVals.a[2]       = mxGetScalar(prhs[10]);
        apVals.a[3]       = mxGetScalar(prhs[11]);
        apVals.a[4]       = mxGetScalar(prhs[12]);
        apVals.a[5]       = mxGetScalar(prhs[13]);
        apVals.a[6]       = mxGetScalar(prhs[14]);
        input.ap_a   = &apVals;
    }

    input.doy    = doy_in;
    input.sec    = sec_in;
    input.alt    = alt_in;
    input.g_lat  = lat_in;
    input.g_long = lon_in;
    input.lst    = lst_in;
    input.f107A  = f107A_in;
    input.f107   = f107_in;
    input.ap     = ApDaily_in;

    gtd7(&input, &flags, &output);

    // Output from MSIS
    plhs[0] = mxCreateDoubleScalar(output.d[0]);
    plhs[1] = mxCreateDoubleScalar(output.d[1]);
    plhs[2] = mxCreateDoubleScalar(output.d[2]);
    plhs[3] = mxCreateDoubleScalar(output.d[3]);
    plhs[4] = mxCreateDoubleScalar(output.d[4]);
    plhs[5] = mxCreateDoubleScalar(output.d[5]);
    plhs[6] = mxCreateDoubleScalar(output.d[6]);
    plhs[7] = mxCreateDoubleScalar(output.d[7]);
    plhs[8] = mxCreateDoubleScalar(output.d[8]);
    plhs[9] = mxCreateDoubleScalar(output.t[0]);
    plhs[10] = mxCreateDoubleScalar(output.t[1]);
}

