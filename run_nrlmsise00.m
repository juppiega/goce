% This runs the NRLMSISE00 model as a height profile for a given time &
% place.
% 
% INPUT:
% h             height vector [km]
% t_datenum     time in datenum format
% glat          geographic latitude [deg]
% glon          geographic longitude [deg]
% F107A         81 day average of F10.7 solar flux
% F107          daily F10.7 solar flux for previous day
% ap            Daily magnetic index
%
%
% OUTPUT: Columns of the out profile are
% 
% 1 height (equal to the input but handy to have)
% 2 HE NUMBER DENSITY(CM-3)
% 3 O NUMBER DENSITY(CM-3)
% 4 N2 NUMBER DENSITY(CM-3)
% 5 O2 NUMBER DENSITY(CM-3)
% 6 AR NUMBER DENSITY(CM-3)                       
% 7 TOTAL MASS DENSITY(GM/CM3) [includes d[8] in td7d]
% 8 H NUMBER DENSITY(CM-3)
% 9 N NUMBER DENSITY(CM-3)
% 10 Anomalous oxygen NUMBER DENSITY(CM-3)
% 11 EXOSPHERIC TEMPERATURE
% 12 TEMPERATURE
%
% NEEDS: 
% "nrlmsise-00.c" and "nrlmsise-00_data.c" from http://www.brodo.de/space/nrlmsise/ 
% and nrlmsise_mex.c written by Monika.
%
% COMPILING: mex CFLAGS="\$CFLAGS -std=c99" COPTIMFLAGS="-O3" nrlmsise_mex.c nrlmsise-00.c nrlmsise-00_data.c
% 
% 
% Monika Andersson & Antti Kero 2013

function out = run_nrlmsise00(h,t_datenum,glat,glon,F107A,F107,ApDaily,apNow,ap3h,ap6h,ap9h,apAver12To33h,apAver36To57h)


    % calculate the day of the year
    year=str2double(datestr(t_datenum,'yyyy'));
    doy=fix(t_datenum)-datenum(year,1,1)+1;
   
    % second of the day
    seconds=(t_datenum-fix(t_datenum))*24*60*60;
    
    % calculate the solar apparent time in hours (roughly)
    lst=seconds/(60*60)+12*glon/180;  % this is just the "local time", in the longitudinal sense, recommended in the model's description.
    
    out=nan(length(h),12); % initialise the output
    out(:,1)=h; % let's put the height to the first column
 
    for hind=1:length(h)
          [out(hind,2),out(hind,3),out(hind,4),out(hind,5),out(hind,6),out(hind,7),out(hind,8),out(hind,9),out(hind,10),out(hind,11),out(hind,12)]...
              =nrlmsise_mex(doy,seconds,h(hind),glat,glon,lst,F107A,F107,ApDaily,apNow,ap3h,ap6h,ap9h,apAver12To33h,apAver36To57h);
    end

return
