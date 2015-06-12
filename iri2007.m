function out = iri2007(time, latitude, longitude, height, utc, coord, ...
    curldir, Rz12, IG12, tec_hmax, ne_top, fpeak, f2storm, bottom, ...
    f1_prob, Ne_Dreg, Te_top, ioncomp, NmF2_foF2, hmF2_M3000F2)

% IRI2007 International Reference Ionosphere 2007 output parameters.
% 
% Usage: OUT = IRI2007(TIME, LATITUDE, LONGITUDE, HEIGHT, UTC, COORD,
%                  CURLDIR, Rz12, IG12, TEC_HMAX, Ne_TOP, FPEAK, F2STORM,
%                  BOTTOM, F1_PROB, Ne_Dreg, Te_TOP, IONCOMP, NmF2_foF2,
%                  hmF2_M3000F2)
% 
% Computes the International Reference Ionosphere (IRI) 2007, which is an
% internationally recognized model for various ionospheric properties. The
% position and time inputs can be scalars or arrays; when they are arrays,
% they should all have the same number of elements. The output is a matrix
% with 40 columns and the same number of rows as elements in the position
% and time inputs (possibly just 1). The ith row has the IRI output for the
% ith position/time input.
% 
% The function makes the computation by querying the online interface at
% http://omniweb.gsfc.nasa.gov/vitmo/iri_vitmo.html (hence internet access
% is required), which makes it pretty slow, especially when either more
% than one of the position and time inputs are made to vary or if the
% position and time inputs are not spaced linearly. If more than one input
% is varied, the function can be sped up by using a for loop and holding
% the smaller arrays constant, assuming the largest array is spaced
% linearly. This will result in fewer calls to the website since the
% website allows for a linear sweep in one variable. See the script
% iritest.m for an example of this. There is one exception: The online
% interface has odd behavior that does not allow for sweeps in longitude
% for any altitude except the default (100 km), so longitude sweeps will be
% computed one at a time unless ALTITUDE is 100.
% 
% The query is made using the command curl in an operating system terminal.
% This is built-in to Unix but not Windows. curl for Windows can be
% downloaded from http://curl.haxx.se/download.html. The directory where
% the curl.exe file can be found should be passed into CURLDIR for Windows
% computers. CURLDIR defaults to the same directory as this function.
% 
% A value for -1 means the output is invalid for the given input. A value
% of -1 for TEC means the input TEC_HMAX was not set. All string inputs are
% case insensitive.
% 
% Inputs:
%   
%   -TIME: Times to compute IRI model either in MATLAB serial date number
%   format or a string that can be converted into MATLAB serial date number
%   format using DATENUM with no format specified (see documentation of
%   DATENUM for more information). Whether the times are local or UTC are
%   determined by the input UTC. Valid range is from year 1958 to year 2015
%   currently (optional, default is January 1, 2000 at 01:30).
%   
%   -LATITUDE: Latitude in degrees to compute IRI model. Whether this is
%   geodetic, geocentric, or geomagnetic latitude is determined by the
%   input COORD. Valid range is -90 degrees to 90 degrees (optional,
%   default is 50 degrees).
%   
%   -LONGITUDE: Longitude in degrees to compute IRI model. Whether this is
%   geodetic, geocentric, or geomagnetic longitude is determined by the
%   input COORD. Valid range is 0 degrees to 360 degrees (optional, default
%   is 40 degrees).
%   
%   -HEIGHT: For geodetic or geomagnetic coordiates, the height in km above
%   the Earth's surface. For geocentric coordiates, the radius in km from
%   the center of the Earth. Valid range for altitude is 60 km to 2000 km,
%   although the recommended upper limit is 1500 km (optional, default when
%   all other inputs are scalars is to sweep from 100:50:2000 km and when
%   any other input is an array is 100 km).
%   
%   -UTC: Set to true (or 'UTC', 'U') if the times in TIME are in
%   Coordinated Universal Time (UTC) or false (or 'Local', 'LT', 'L') if
%   the times in TIME are in local time (optional, default is true).
%   
%   -COORD: String specifying the coordinate system to use. Can be
%   geodetic ('geodetic', 'geod', or 'gd'), geomagnetic ('geomagnetic',
%   'geom', or 'gm'), or geocentric ('geocentric', 'geom', or 'gm')
%   (optional, default is geodetic).
%   
%   -CURLDIR: Directory where curl.exe can be found (optional, only
%   necessary for Windows computers, and default for those is the same
%   directory that this function is located).
%   
%   -Rz12: 13-month running mean of sunspot number index. This index is
%   typically calculated automatically from the time you specify, but you
%   can enter your own index here and not rely on the internal index file.
%   Valid range is 0 to 400.
%   
%   -IG12: Index based on foF2 measurements from a dozen ionosondes
%   correlated with the CCIR foF2 maps [1].  This index is typically
%   calculated automatically from the time you specify, but you can enter
%   your own index here and not rely on the internal index file. Valid
%   range is -50 to 400.
%   
%   -TEC_HMAX: Electron content upper boundary in km. The electron content
%   is the vertical electron content (vTEC) calculated using numerical
%   integration from a lower to an upper boundary. You must enter a value
%   for the upper boundary if you want to obtain electron content values.
%   The value should not be higher than 2000 km because this is the upper
%   limit of the IRI validity range for electron density profiles. The
%   lower boundary is set to the starting height of the IRI validity range
%   (60 km during the day and 80 km during the night).
% 
%   -Ne_TOP: Ne Topside. There are three options here:
%     1. IRI-2001: The IRI-2001 model that was based primarily on Alouette
%     1 topside sounder data with some AE-C, AEROS, and DE-2 in situ data
%     [2]. Selected by inputting one of 3, 'IRI01', 'IRI2001', or '2001'.
%     2. IRI01-corr: A correction of the 2001 model with the help of
%     Aloutte 2, ISIS 1 and 2 topside sounder data [3]. Selected by
%     inputting one of 2, 'IRI01-corr', 'IRI2001-corr', 'corr', or 'c'.
%     3. NeQuick: The model developed by Radicella and his group at ICTP
%     using Intercosmos 19 topside sounder data in addition to the ISIS 1
%     and 2 data [4]. Selected by inputing one of 1, 'NeQuick', 'Quick',
%     'Qu', 'Q', or 'Ne'.
%   The default option is NeQuick.
% 
%   -FPEAK: F peak model. There are two options here:
%     1. CCIR: The CCIR model for the F peak plasma frequency foF2 was
%     developed by Jones and Gallet [5] using data from the worldwide
%     network of ionosondes. It is the model recommended by the Comite
%     Consultatif International des Radiocommunications (CCIR) of the
%     International Telecommunications Union [6]. Selected by inputting one
%     of 2, 'CCIR', or 'C'.
%     2. URSI: The URSI model was developed by Rush et al. [7] using a
%     physical model to obtain screen points over regions not covered by
%     ionosondes instead of the extrapolation along magnetic field lines
%     employed by Jones and Gallet [5]. Selected by inputting one of 1,
%     'URSI', or 'U'.
%   The default option is URSI. The CCIR model is recommended over the
%   continents and the URSI model over the oceans. NOTE: Changes here will
%   also affect the height of the F peak hmF2 because the hmF2 model
%   depends on the ratio foF2/foE where foE is the plasma frequency at the
%   E peak.
%   
%   -F2STORM: foF2 Storm model on (true or 'on') or off (false or 'off')
%   (optional, default is on). The F peak storm model describes the average
%   storm behavior in terms of the ratio foF2_storm/foF2_quiet based on the
%   ap history over the preceding 33 hours [8]. A large volume of ionosonde
%   data for storms during the 1980-1990 time period were used to describe
%   the most coherent and repeatable features of the ionospheric storm
%   response.
% 
%   -BOTTOM: Bottomside thickness. The bottomside thickness is the height
%   difference between the F peak height hmF2 and the height where the
%   electron density profile has dropped down to half the F peak value
%   (NmF2/2). There are three options here:
%     1. Bilitza-2000: This (table-)option is based on incoherent scatter
%     data [9]. Selected by inputting one of 1, 'B0 Table', 'B0Table',
%     'B0', or 'B'.
%     2. Gulyaeva: This option is based on ionosonde data mostly from mid-
%     latitudes [10]. Selected by inputting one of 2, 'Gulyaeva', 'Guly',
%     or 'G'.
%   The default option is Bilitza-2000.
% 
%   -F1_PROB: This parameter describes the occurrence probability of an F1
%   layer. There are three options here:
%     1. IRI-95: This option uses the ITU-recommended Ducharme et al. model
%     [11] that applies a simple cutoff solar zenith angle, so probability
%     is either 0 or 1. Selected by inputting one of 3, 'IRI-95', 'IRI95',
%     'IRI', or '95'.
%     2. Scotto-1997 no L: This option uses the model developed by Scotto
%     et al. [12] using only ionograms with clear F1 layer presence and
%     excluding the more uncertain L condition cases. Ionograms often
%     exhibit a Fl ledge rather than a fully developed cusp, primarily
%     during the time period just before the Fl layer disappears. These
%     cases are described as L condition according to the URSI standard
%     nomenclature. Selected by inputting one of 1, 'Scotto-1997 no L', 
%     'no L', or 'nL'.
%     3. Scotto-1997 with L: Here L-condition cases were included. Selected
%     by inputting one of 2, 'Scotto-1997 with L', 'with L', or 'wL'.
%   The default option is Scotto-1997 no L.
% 
%   -Ne_Dreg: D-region electron density. There are two options here:
%     1. FT-2001: Model based on Friedrich's compilation of reliable rocket
%     data [13]. Selected by inputting one of 2, 'FPT-2000', 'FPT', or
%     '2000'.
%     2. IRI-95: Model based on a much smaller selection of typical rocket
%     profiles [14]. Selected by inputting one of 1, 'IRI-95', 'IRI95',
%     'IRI', or '95'.
%   The default option is IRI-95.
% 
%   -Te_TOP: Topside electron temperature. There are two options here:
%     1. IRI-95 is using the global models at fixed altitudes developed by
%     Brace and Theis [15] based on their ISIS and AE-C electron
%     temperature measurements. The model is described in [16] and also in
%     the IRI-90 report that is available as PDF document from the
%     references section of the IRI homepage. Selected by inputting one of
%     2, 'IRI-95', 'IRI95', 'IRI', or '95'.
%     2. TTSA-2000 model is newer and is based primarily on Intercosmos 19,
%     24, and 25 data [17, 18]. Selected by inputting one of 1,
%     'TTSA-2000', 'TTSA', or '2000'.
%   The default option is TTSA-2000.
% 
%   -IONCOMP: Ion composition. There are two options here:
%     1. DS95/TTS05: Uses the Danilov and Smirnova model [19] based on
%     their compilation of rocket data in the region below the F-peak and
%     the T?ísková, Truhlík, and Šmilauer model [20] based on satellite
%     data in the region above the peak. Selected by inputting one of 1,
%     'DS95/TTS05', 'DS95', or '95'.
%     2. DS78/DY85: The pre-2000 model used in IRI based on the work of
%     Danilov and Semenov [21] and of Danilov and Yaichnikov [22] with
%     their selection of high-alittutde rocket data. Selected by inputting
%     one of 2, 'DS78/DY95', or 'DS78', or '78'.
%   The default option is DS95/TTS05.
%   
%   -NmF2_foF2: A measured value to update the IRI profile to actual
%   conditions. Only valid when varying HEIGHT linearly. If this number is
%   between 10^3 and 10^8, it represents the F2 peak density (NmF2) in
%   cm^-3. If instead it is between 2 and 14, it represents the F2 plasma
%   frequency foF2 in MHz.
%   
%   -hmF2_M3000F2: A measured value to update the IRI profile to actual
%   conditions. Only valid when varying HEIGHT linearly. If this number is
%   between 100 and 1000, it represents the F2 peak height (hmF2) in km. If
%   instead it is between 1.5 and 4, it represents the propagation factor
%   M(3000)F2.
% 
% Outputs:
%   
%   -OUT: Array with the same number of rows as elements in the position
%   and time inputs and the following data in each column:
%     1. Electron density (Ne) in m^-3.
%     2. Ratio of Ne to the F2 peak density (NmF2).
%     3. Neutral temperature (Tn) in K.
%     4. Ion temperature (Ti) in K.
%     5. Electron temperature (Te) in K.
%     6. Atomic oxygen ions (O+) percentage.
%     7. Atomic hydrogen ions (H+) percentage.
%     8. Atomic helium ions (He+) percentage.
%     9. Molecular oxygen ions (02+) percentage.
%     10. Nitric oxide ions (NO+) percentage.
%     11. Cluster ions percentage.
%     12. Atomic nitrogen ions (N+) percentage.
%     13. Total electron content (TEC) in 10^16 m^-2.
%     14. TEC top percentage.
%     15. Height of the F2 peak (hmF2) in km.
%     16. Height of the F1 peak (hmF1) in km.
%     17. Height of the E peak (hmE) in km.
%     18. Height of the D peak (hmD) in km.
%     19. Density of the F2 peak (NmF2) in m^-3.
%     20. Density of the F1 peak (NmF1) in m^-3.
%     21. Density of the E peak (NmE) in m^-3.
%     22. Density of the D peak (NmD) in m^-3.
%     23. Propagation factor M(3000)F2.
%     24. Bottomside thickness (B0) in km.
%     25. Bottomside shape (B1).
%     26. E-valley width in km.
%     27. E-valley depth (Nmin/NmE) in km.
%     28. F2 plasma frequency (foF2) in MHz.
%     29. F1 plasma frequency (foF1) in MHz.
%     30. E plasma frequency (foE) in MHz.
%     31. D plasma frequency (foD) in MHz.
%     32. Equatorial vertical ion drift in m/s.
%     33. Ratio of foF2 storm to foF2 quiet.
%     34. F1 probability 1979.
%     35. F1 probability 1995 including L-condition.
%     36. F1 probability 1995.
%     37. Spread F probability.
%     38. 12-month running mean of sunspot number Rz12 used by the model.
%     39. Ionospheric index IG12 used by the model.
%     40. Daily solar radio flux F107D used by the model.
%   A value of -1 for any output indicates that the parameter is not
%   available for the specified range. TEC = -1 means you have not entered
%   an upper boundary height for TEC_HMAX.
% 
% References:
% 
%   [1] R. Y. Liu, P. A. Smith, and J. W. King, “A New Solar Index Which
%   Leads to Improved foF2 Predictions Using the CCIR Atlas,"
%   _Telecommunication_Journal_, vol. 50, no. 8, pp. 408-414, 1983.
%   [2] D. Bilitza, B. W. Reinisch, S. M. Radicella, S. Pulinets, T.
%   Gulyaeva, and L. Triskova, "Improvements of the International Reference
%   Ionosphere model for the topside electron density profile," _Radio_
%   _Sci._, vol. 41, p. RS5S15, 2006, doi:10.1029/2005RS003370.
%   [3] D. Bilitza, "A correction for the IRI topside electron density
%   model based on Alouette/ISIS topside sounder data," _Adv._Space_Res._,
%   vol. 33, no. 6, 838-843, 2004.
%   [4] Coisson P., S. Radicella, R. Leitinger, and B. Nava, "Topside
%   electron density in IRI and NeQuick," _Adv._Space_Res._, vol. 37, no.
%   5, pp. 937-942, 2006.
%   [5] W. B. Jones and R. M. Gallet, "The Representation of Diurnal and
%   Geographic Variations of Ionospheric Data by Numerical Methods,"
%   _Telecomm._J._, vol. 29, no. 129, 1962, and vol. 32, no. 18, 1965. 
%   [6] Report 340-4, ITU, Geneva, 1967.
%   [7] C. Rush, M. Fox. D. Bilitza, K. Davies, L. McNamara, F. Stewart,
%   and M. PoKempner, "Ionospheric Mapping-An Update of foF2 Coefficients,"
%   _Telecomm._J._, vol. 56, no. 179, 1989.
%   [8] T. Fuller-Rowell, E. Araujo-Pradere, and M. Codrescu, "An empirical
%   ionospheric storm-time correction model," _Adv._Space_Res._, vol. 25,
%   no. 1, pp. 139-148, 2000.
%   [9] D. Bilitza, S. M. Radicella, B. W. Reinisch, J. O. Adeniyi, M. E.
%   Mosert Gonzalez, S. R. Zhang, O. Obrou, "New B0 and B1 models for IRI,"
%   _Adv._Space_Res._, vol. 25, no. 1, pp. 89-95, 2000,
%   doi:10.1016/S0273-1177(99)00902-3.
%   [10] T. L. Gulyaeva, "Progress in ionospheric informatics based on
%   electron-density profile analysis of ionograms," _Adv._Space_Res._,
%   vol. 7, no. 6, 1987, doi:10.1016/0273-1177(87)90269-9.
%   [11] E. D. DuCharme, L. E. Petrie, and R. Eyfrig, "A method for
%   predicting the F1 layer critical frequency based on the Zurich smoothed
%   sunspot number," _Radio_Sci._, vol. 8, no. 10, pp. 837-839, 1973,
%   doi:10.1029/RS008i010p00837.
%   [12] C. Scotto, S. M. Radicella, and B. Zolesi, "An improved
%   probability function to predict the F1 layer occurrence and L
%   condition," _Radio_Sci._, vol. 33, no. 6, pp. 1763-1765, 1998,
%   doi:10.1029/98RS02637.
%   [13] M. Friedrich and K. M. Torkar, "FIRI: A semiempirical model of the
%   lower ionosphere," _J._Geophys._Res._, vol 106, no. A10, 2001,
%   doi:10.1029/2001JA900070.
%   [14] E. Mechtly and D. Bilitza, "Models of D-region electron
%   concentration," Rep. IPW-WB1, Inst. fur phys. Weltraumforsch.,
%   Freiburg, Germany, 1974.
%   [15] L. H. Brace and R. F. Theis, "Global empirical models of
%   ionospheric electron temperature in the upper F region and plasmasphere
%   based on in situ measurements from the Atmospheric Explorer-C, ISIS-1
%   and ISIS-2 satellites," _J._Atmos._Terr._Phys._, vol. 43, pp.
%   1317–1343, 1981.
%   [16] D. Bilitza, L. H. Brace, and R. F. Theis, "Modelling of
%   ionospheric temperature profiles," _Adv._Space_Res._, vol. 5, no. 7,
%   pp. 53-58, 1985, doi:10.1016/0273-1177(85)90356-4.
%   [17] V. Truhlík, L. T?ísková, J. Šmilauer, and V. V. Afonin, "Global
%   empirical model of electron temperature in the outer ionosphere for
%   period of high solar activity based on data of three Intercosmos
%   satellites," _Adv._Space_Res._, vol. 25, no. 1, pp. 163-169, 2000 
%   doi:10.1016/S0273-1177(99)00914-X.
%   [18] V. Truhlík, L. T?ísková, J. Šmilauer, and V. V. Afonin, "Improved
%   electron temperature model and comparison with satellite data," 
%   _Adv._Space_Res._, vol. 27, no. 1, pp. 101-109, 2001
%   doi:10.1016/S0273-1177(00)00144-7.
%   [19] A. D. Danilov, and N. V. Smirnova, "Improving the 75 to 300 km ion
%   composition model of the IRI," _Adv._Space_Res._, vol. 15, no. 2, pp.
%   171-177, 1995, doi:10.1016/S0273-1177(99)80044-1.
%   [20] L. T?ísková, V. Truhlík,, and J. Šmilauer, "An empirical model of
%   ion composition in the outer ionosphere," _Adv._Space_Res._, vol. 31,
%   no. 3, pp. 653-663, 2003, doi:10.1016/S0273-1177(03)00040-1.
%   [21] A. D. Danilov, V. K. Semenov, "Relative ion composition model at
%   midlatitudes," vol. 40, nos. 10-11, pp. 1093-1102, 1978,
%   doi:10.1016/0021-9169(78)90057-0.
%   [22] A. D. Danilov, and A. P. Yaichnikov, "A new model of the ion
%   composition model at 75 to 1000 km for IRI," _Adv._Space_Res._, 
%   vol. 5, no. 7, pp. 75–79, 1985, doi:10.1016/0273-1177(85)90360-6.
% 
% See also: IRI2012, IRITEST, MSIS, IGRF, DATENUM.

% Directory of this function.
fpath = fileparts(mfilename('fullpath'));

% Default behavior.
if nargin < 1 || isempty(time)
    time = datenum([2000 1 1 1 30 0]);
end
if nargin < 2 || isempty(latitude)
    latitude = 50;
end
if nargin < 3 || isempty(longitude)
    longitude = 40;
end
if nargin < 4 || isempty(height)
    if (ischar(time) || numel(time) == 1) && ...
            numel(latitude) == 1 && numel(longitude) == 1
        height = 100:50:2000;
    else
        height = 100;
    end
end
if nargin < 5 || isempty(utc)
    utc = true;
elseif ischar(utc)
    switch lower(utc)
        case {'utc', 'u'}
            utc = true;
        case {'local', 'lt', 'l'}
            utc = false;
        otherwise
            error('iri:badUTC', ['Unrecognized command ' utc '. Valid ' ...
                'options are ''utc'', ''u'', ''local'', ''lt'', or ' ...
                '''l''.']);
    end
end
if nargin < 6 || isempty(coord)
    coord = 'geod';
end
if isunix || ismac  % Curl is built-in to operating system.
    curldir = [];
elseif nargin < 7 || isempty(curldir)
    curldir = fpath;
end
if nargin < 8
    Rz12 = [];
end
if nargin < 9
    IG12 = [];
end
if nargin < 10
    tec_hmax = [];
end
if nargin < 11 || isempty(ne_top)
    ne_top = 'NeQuick';
end
if nargin < 12 || isempty(fpeak)
    fpeak = 'URSI';
end
if nargin < 13 || isempty(f2storm)
    f2storm = true;
elseif ischar(f2storm)
    if strcmpi(f2storm, 'on')
        f2storm = true;
    elseif strcmpi(f2storm, 'off')
        f2storm = false;
    else
        error('iri:invalidf2storm', ...
            'Unrecognized input ''%s'' for F2STORM.', f2storm);
    end
end
if nargin < 14 || isempty(bottom)
    bottom = 'B0 Table';
end
if nargin < 15 || isempty(f1_prob)
    f1_prob = 'Scotto-1997 no L';
end
if nargin < 16 || isempty(Ne_Dreg)
    Ne_Dreg = 'IRI-95';
end
if nargin < 17 || isempty(Te_top)
    Te_top = 'TTSA-2000';
end
if nargin < 18 || isempty(ioncomp)
    ioncomp = 'DS95/TTS05';
end
if nargin < 19 || isempty(NmF2_foF2)
    NmF2_foF2 = 0;
end
if nargin < 20 || isempty(hmF2_M3000F2)
    hmF2_M3000F2 = 0;
end

% Convert time to a datenumber if it is a string.
if ischar(time)
    time = datenum(time);
end

% Convert the coordinates.
switch lower(coord)
    case {'geodetic', 'geod', 'gd'}
        geo_flag = '0.';
    case {'geocentric', 'geoc', 'gc'}
        % Convert coordinates to geodetic. The function ecef2geod assumes
        % meters, but we want km here.
        [x, y, z] = sph2cart(longitude*pi/180, latitude*pi/180, ...
            height*1e3);
        [latitude, longitude, height] = ecef2geod(x, y, z);
        height = height/1e3;
        geo_flag = '0.';
    case {'geomagnetic', 'geom', 'gm'}
        geo_flag = '1.';
    otherwise
        error('iri:coordCommandUnknown', ['Command ' coord ' unknown. ' ...
            'Valid options are ''geodetic'', ''geocentric'', and ' ...
            '''geomagnetic''.']);
end

% Error checking and input conversion.
if any(latitude < -90) || any(latitude > 90)
    error('iri:invalidLatitude', ['Input LATITUDE must be between ' ...
        '-90 degrees and 90 degrees.']);
end
latitude = latitude(:).';
longitude = mod(longitude(:).', 360);
if any(height < 0) || any(height > 2e3)
    error('iri:invalidHeight', ...
        'Input HEIGHT must be between 0 km and 2000 km.');
end
height = height(:).';
if isempty(Rz12) || (Rz12 >= 0 && Rz12 <= 400)
    Rz12 = sprintf('&sun_n=%#g', Rz12);
else
    error('iri:invalidRz12', ...
        'Input RZ12 is %g but must be between 0 and 400.', Rz12);
end
if isempty(IG12) || (IG12 >= -50 && IG12 <= 400)
    IG12 = sprintf('&ion_n=%#g', IG12);
else
    error('iri:invalidIG12', ...
        'Input IG12 is %g must be between -50 and 400.', IG12);
end
if isempty(tec_hmax) || (tec_hmax >= 50 && tec_hmax <= 2000)
    tec_hmax = sprintf('&htec_max=%#g', tec_hmax);
else
    error('iri:invalidTec_hmax', ['Input TEC_HMAX is %g km but must ' ...
        'be between 50 km and 2000 km.'], tec_hmax);
end
switch lower(ne_top)
    case {1, 'nequick', 'quick', 'qu', 'q', 'ne'}
        ne_top = '0.';
    case {2, 'iri01-corr', 'iri2001-corr', 'corr', 'c'}
        ne_top = '1.';
    case {3, 'iri01', 'iri2001', '2001'}
        ne_top = '2.';
    otherwise
        if ~ischar(ne_top)
            ne_top = num2str(ne_top);
        else
            ne_top = ['' ne_top ''];
        end
        error('iri:invalidNe_top', ['Unrecognized command ' ne_top ...
            ' for input NE_TOP. Valid options are NeQuick (1, ' ...
            '''NeQuick'', ''Quick'', ''Qu'', ''Q'', or ''Ne''), ' ...
            'IRI01-corr (2, ''IRI01-corr'', ''IRI2001-corr'', ' ...
            '''corr'', or ''c'') and IRI2001 (3, ''IRI01'', ' ...
            '''IRI2001'', or ''2001'').']);
end
switch lower(fpeak)
    case {1, 'ursi', 'u'}
        fpeak = '0.';
    case {2, 'ccir', 'c'}
        fpeak = '1.';
    otherwise
        if ~ischar(fpeak)
            fpeak = num2str(fpeak);
        else
            fpeak = ['' fpeak ''];
        end
        error('iri:invalidFpeak', ['Unrecognized command ' fpeak ...
            ' for input FPEAK. Valid options are URSI (1, ''URSI'', ' ...
            ' or ''U'') and CCIR (2, ''CCIR'', or ''C'').']);
end
if f2storm
    f2storm = '0.';
else
    f2storm = '1.';
end
switch lower(bottom)
    case {1, 'b0 table', 'b0table', 'b0', 'b'}
        bottom = '0.';
    case {2, 'gulyaeva', 'guly', 'g'}
        bottom = '1.';
    otherwise
        if ~ischar(bottom)
            bottom = num2str(bottom);
        else
            bottom = ['' bottom ''];
        end
        error('iri:invalidBottom', ['Unrecognized command ' bottom ...
            ' for input BOTTOM. Valid options are Bilitza-2000 (1, ' ...
            '''B0 Table'', ''B0Table'', ''B0'', or ''B''), and ' ...
            'Gulyaeva (2, ''Gulyaeva'', ''Guly'', or ''G'').']);
end
switch lower(f1_prob)
    case {1, 'scotto-1997 no l', 'no l', 'nl'}
        f1_prob = '0.';
    case {2, 'scotto-1997 with l', 'with l', 'wl'}
        f1_prob = '1.';
    case {3, 'iri-95', 'iri95', 'iri', '95'}
        f1_prob = '2.';
    otherwise
        if ~ischar(f1_prob)
            f1_prob = num2str(f1_prob);
        else
            f1_prob = ['' f1_prob ''];
        end
        error('iri:invalidF1prob', ['Unrecognized command ' f1_prob ...
            ' for input F1PROB. Valid options are Scotto-1997 no L ' ...
            '(1, ''Scotto-1997 no L'', ''no L'', or ''nL'') ' ...
            'Scotto-1997 with L (2, ''Scotto-1997 with L'', ' ...
            '''with L'', or ''wL''), and IRI-95 (3, ''IRI-95'', ' ...
            '''IRI95'', ''IRI'', or ''95'').']);
end
switch lower(Ne_Dreg)
    case {1, 'iri-95', 'iri95', 'iri', '95'}
        Ne_Dreg = '0.';
    case {2, 'fpt-2000', 'fpt', '2000'}
        Ne_Dreg = '1.';
    otherwise
        if ~ischar(Ne_Dreg)
            Ne_Dreg = num2str(Ne_Dreg);
        else
            Ne_Dreg = ['' Ne_Dreg ''];
        end
        error('iri:invalidD_reg', ['Unrecognized command ' Ne_Dreg ...
            ' for input D_REG. Valid options are IRI-95 (1, ' ...
            '''IRI-95'', ''IRI95'', ''IRI'', or ''95'') and ' ...
            'FPT-2000 (2, ''FPT-2000'', ''FPT'', or ''2000'').']);
end
switch lower(Te_top)
    case {1, 'ttsa-2000', 'ttsa', '2000'}
        Te_top = '0.';
    case {2, 'iri-95', 'iri95', 'iri', '95'}
        Te_top = '1.';
    otherwise
        if ~ischar(Te_top)
            Te_top = num2str(Te_top);
        else
            Te_top = ['' Te_top ''];
        end
        error('iri:invalidTe_top', ['Unrecognized command ' Te_top ...
            ' for input TE_TOP. Valid options are TTSA-2000 (1, ' ...
            '''TTSA-2000'', ''TTSA'', or ''2000'') and IRI-95 (2, ' ...
            '''IRI-95'', ''IRI95'', ''IRI'', or ''95'').']);
end
switch lower(ioncomp)
    case {1, 'ds95/tts05', 'ds95', '95'}
        ioncomp = '0.';
    case {2, 'ds78/dy85', 'ds78', '78'}
        ioncomp = '1.';
    otherwise
        if ~ischar(ioncomp)
            ioncomp = num2str(ioncomp);
        else
            ioncomp = ['' ioncomp ''];
        end
        error('iri:invalidIoncomp', ['Unrecognized command ' ioncomp ...
            ' for input IONCOMP. Valid options are DS95/TTS05 (1, ' ...
            '''DS95/TTS05'', ''DS95'', or ''95'') and DS78/DY85 (2, ' ...
            '''DS78/DY85'', ''DS78'', or ''78'').']);
end
if ~(NmF2_foF2 == 0 || (NmF2_foF2 >= 2 && NmF2_foF2 <= 14) || ...
        (NmF2_foF2 >= 1e3 && NmF2_foF2 <= 1e8))
    error('iri:invalidNmF2_foF2', ['Input NmF2_foF2 is %g but must be ' ...
        'one of 0, between 2 and 14 (MHz), or between 10^3 and 10^8 ' ...
        '(cm^-3).'], ...
        NmF2_foF2);
end
if ~(hmF2_M3000F2 == 0 || ...
        (hmF2_M3000F2 >= 100 && hmF2_M3000F2 <= 1000) || ...
        (hmF2_M3000F2 >= 1.5 && hmF2_M3000F2 <= 4))
    error('iri:invalidhmF2_M3000F2', ['Input hmF2_M3000F2 is %g but ' ...
        'must be one of 0, between either 1.5 and 4, or between 100 ' ...
        'and 1000 (km).'], hmF2_M3000F2);
end

% Get the file name to create for the curl command. To allow function to
% run in a parfor loop, append the current task and a random number to the
% usual file name.
if license('test', 'Distrib_Computing_Toolbox')
    t = getCurrentTask();
    if isempty(t)
        fname = 'temp';
    else
        fname = sprintf('temp%i', t.ID);
    end
else
    fname = 'temp';
end
fname = fullfile(fpath, sprintf('%s_%i.html', fname, randi(2^53-1, 1)));

% End of the string to be input into the system function.
endcmd = [Rz12 IG12 tec_hmax ...
    '&ne_top=' ne_top ...
    '&imap=' fpeak ...
    '&ffof2=' f2storm ...
    '&ib0=' bottom ...
    '&probab=' f1_prob ...
    '&dreg=' Ne_Dreg ...
    '&tset=' Te_top ...
    '&icomp=' ioncomp ...
    sprintf('&nmf2=%#g', NmF2_foF2) ...
    sprintf('&f2=%#g', hmF2_M3000F2) ...
    sprintf('&vars=%i', [6, 11:50]) ...
    '" http://omniweb.gsfc.nasa.gov/cgi/vitmo/vitmo_model.cgi > "' ...
    ...'" http://ccmc.gsfc.nasa.gov/cgi-bin/modelweb/models/vitmo_model.cgi > "' ...
    fname '"'];

% Get year, month, day, hour for IRI.
[year, month, day, hour, minute, second] = datevec(time);
hour = hour + minute/60 + second/3600;
dayyear = dayofyear(year, month, day);

% Get the beginning of the string to be input into the system function.
if isunix || ismac % Curl is built-in to operating system.
    initialcmd = 'curl -d "model=iri';
else % Use the input curldir (possibly the default for that input).
    initialcmd = ['"' fullfile(curldir, 'curl"') ' -d "model=iri'];
end

% If one of the inputs can be swept keeping the others constant, run a
% sweep on the model.
A = [numel(unique(height)), numel(unique(latitude)), ...
    numel(unique(longitude)), numel(unique(year)), ...
    numel(unique(dayyear)), numel(unique(hour))];
if all(A == [1 1 1 1 1 1]) % Just one unique input.
    profile = 0;
elseif all(A == [A(1) 1 1 1 1 1]) % Sweep altitude.
    [sweep, tmp, sortind] = unique(height);
    if all(abs(diff(diff(sweep))) < 1e-12)
        profile = 1;
        height = sweep(1);
        sweepmax = 2000;
    else % Altitude is unsweepable because steps are not linear.
        profile = 9;
    end
elseif all(A == [1 A(2) 1 1 1 1]) % Sweep latitude.
    [sweep, tmp, sortind] = unique(latitude);
    if all(abs(diff(diff(sweep))) < 1e-12)
        profile = 2;
        latitude = sweep(1);
        sweepmax = 90;
    else % Latitude is unsweepable because steps are not linear.
        profile = 9;
    end
elseif all(A == [1 1 A(3) 1 1 1]) % Sweep longitude.
    [sweep, tmp, sortind] = unique(longitude);
    if all(abs(diff(diff(sweep))) < 1e-12)
        profile = 3;
        longitude = sweep(3);
        sweepmax = 360;
    else % Longitude is unsweepable because steps are not linear.
        profile = 9;
    end
elseif all(A == [1 1 1 A(4) 1 1]) % Sweep year.
    [sweep, tmp, sortind] = unique(year);
    if all(abs(diff(diff(sweep))) < 1e-12)
        profile = 4;
        year = sweep(1);
        sweepmax = 2016;
    else % Year is unsweepable because steps are not linear.
        profile = 9;
    end
elseif all(A == [1 1 1 1 A(5) 1]) % Sweep day of year or month.
    [sweep, tmp, sortind] = unique(dayyear);
    if all(abs(diff(diff(sweep))) < 1e-12)
        profile = 7;
        dayyear = sweep(1);
        sweepmax = 365 + ( (~mod(year(1), 4)...
             & mod(year(1), 100)) | (~mod(year(1), 400)) );
    % Day of year steps are not linear, but maybe month steps are.
    elseif numel(day(1)) == 1
        [sweep, tmp, sortind] = unique(month);
        profile = 5;
        month = sweep(1);
        sweepmax = 12;
    else % Day of year is unsweepable because steps are not linear.
        profile = 9;
    end
elseif all(A == [1 1 1 1 1 A(6)]) % Sweep hour in a day.
    [sweep, tmp, sortind] = unique(hour);
    if all(abs(diff(diff(sweep))) < 1e-12)
        profile = 8;
        hour = sweep(1);
        sweepmax = 24;
    else % Hour is unsweepable because steps are not linear.
        profile = 9;
    end
else
    profile = 9;
end

% Call curl depending on the sweep profile.
switch profile
    case 0
        [status, result] = system([initialcmd ...
            sprintf('&year=%i', year) ...
            sprintf('&month=%i', month) ...
            sprintf('&day=%i', day) ...
            sprintf('&time_flag=%i', ~utc) ...
            sprintf('&hour=%#g', hour) ...
            '&geo_flag=' geo_flag ...
            sprintf('&latitude=%#g', latitude) ...
            sprintf('&longitude=%#g', longitude) ...
            sprintf('&height=%#g', height) ...
            '&profile=1' ...
            sprintf('&start=%#g', height) ...
            sprintf('&stop=%#g', height) ...
            '&step=1.' ...
            endcmd]);
        if status == 0
            out = parseresult(fname, 1);
        else
            error('iri:curlError', ['Curl command did not work. ' ...
                'It returned status:%i, ' ...
                regexprep(result, '\\', '\\\\')], status);
        end
    case {1 2 3 4 5 7 8}
        % The online interface will only output up to a sweep length of
        % 500, so split the sweep into increments of 500.
        nsweeps = ceil(numel(sweep) / 500);
        sweepstart = sweep(1 : 500 : end);
        sweepstop = sweep(500 : 500 : end);
        if numel(sweepstart) - 1 == numel(sweepstop)
            sweepstop(end + 1) = sweep(end);
        end
        sweepstep = mode(diff(sweep));
        sweeplen = round((sweepstop - sweepstart) ./ sweepstep + 1);
        prevsweeplen = [0 sweeplen(1:end-1)];
        sweepstop = min([sweep(end) + sweepstep/10, sweepmax]);
        out = zeros(41*sum(sweeplen), 1);
        for index = 1:nsweeps
            [status, result] = system([initialcmd ...
                sprintf('&year=%i', year(1)) ...
                sprintf('&month=%i', month(1)) ...
                sprintf('&day=%i', day(1)) ...
                sprintf('&time_flag=%i', ~utc) ...
                sprintf('&hour=%#g', hour(1)) ...
                '&geo_flag=' geo_flag ...
                sprintf('&latitude=%#g', latitude(1)) ...
                sprintf('&longitude=%#g', longitude(1)) ...
                sprintf('&height=%#g', height(1)) ...
                sprintf('&profile=%i', profile) ...
                sprintf('&start=%#g', sweepstart(index)) ...
                sprintf('&stop=%#g', sweepstop) ...
                sprintf('&step=%#g', sweepstep) ...
                endcmd]);
            if status == 0
                out((1:41*sweeplen(index)) + 41*prevsweeplen(index)*...
                    (index-1)) = parseresult(fname, sweeplen(index));
            else
                error('iri:curlError', ['Curl command did not work. ' ...
                    'It returned status:%i, ' ...
                    regexprep(result, '\\', '\\\\')], status);
            end
        end
    % profile = 9 means there is more than one unique run to make. Turn all
    % the scalars into vectors with the same number of elements as the
    % largest array. If there is an array smaller than the largest array,
    % throw an error.
    case 9
        maxnum = max([numel(height), numel(latitude), ...
            numel(longitude), numel(year), numel(day), numel(month), ...
            numel(hour)]);
        if numel(height) == 1
            height = repmat(height, maxnum, 1);
        elseif numel(height) ~= maxnum && numel(height) > 1
            error('iri:invalidSize', ['Input vectors must all have ' ...
                'the same number of elements.']);
        end
        if numel(latitude) == 1
            latitude = repmat(latitude, maxnum, 1);
        elseif numel(latitude) ~= maxnum && numel(latitude) > 1
            error('iri:invalidSize', ['Input vectors must all have ' ...
                'the same number of elements.']);
        end
        if numel(longitude) == 1
            longitude = repmat(longitude, maxnum, 1);
        elseif numel(longitude) ~= maxnum && numel(longitude) > 1
            error('iri:invalidSize', ['Input vectors must all have ' ...
                'the same number of elements.']);
        end
        if numel(time) == 1
            year = repmat(year, maxnum, 1);
            month = repmat(month, maxnum, 1);
            day = repmat(day, maxnum, 1);
            hour = repmat(hour, maxnum, 1);
        elseif numel(time) ~= maxnum && numel(time) > 1
            error('iri:invalidSize', ['Input vectors must all have ' ...
                'the same number of elements.']);
        end
        out = zeros(41*maxnum, 1);
        for index = 1:maxnum
            [status, result] = system([initialcmd ...
                sprintf('&year=%i', year(index)) ...
                sprintf('&month=%i', month(index)) ...
                sprintf('&day=%i', day(index)) ...
                sprintf('&time_flag=%i', ~utc) ...
                sprintf('&hour=%#g', hour(index)) ...
                '&geo_flag=' geo_flag ...
                sprintf('&latitude=%#g', latitude(index)) ...
                sprintf('&longitude=%#g', longitude(index)) ...
                sprintf('&height=%#g', height(index)) ...
                '&profile=1' ...
                sprintf('&start=%#g', height(index)) ...
                sprintf('&stop=%#g', height(index)) ...
                '&step=1.' ...
                endcmd]);
            if status == 0
                out((1:41) + 41*(index-1)) = parseresult(fname, 1);
            else
                error('iri:curlError', ['Curl command did not work. ' ...
                    'It returned status:%i, ' ...
                    regexprep(result, '\\', '\\\\')], status);
            end
        end
end % End case

% Right now, out has all the data in one long vector, but we want it to be
% an array, so reshape.
out(1:41:end) = []; % Remove height data.
out = reshape(out, 40, []).';

% Resort the output if we did a sweep.
if profile >= 1 && profile <= 8
    out = out(sortind, :);
end


% Return the day of the year.
function dayyear = dayofyear(year, month, day)

previous = cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30]);
dayyear = previous(month) + day + double( ...
    (~mod(year, 4) & mod(year, 100)) | (~mod(year, 400)) & (month > 2) );


% Parse the result file output by the system command.
function data = parseresult(fname, sweeplen)

% Get the data from fname into a string, then delete the file.
fid = fopen(fname);
if fid == -1
    error('iri:parseresult:cannotOpenFile', ['Cannot open %s ' ...
        'file generated by curl command.'], fname);
end
result = fread(fid, '*char').';
f = fclose(fid);
if f == -1
    error('iri:parseresult:cannotCloseFile', ['Cannot close %s ' ...
        'file generated by curl command.'], fname);
end
delete(fname);

% Data starts at line 71. If there isn't a line 71, there was an error.
% Also, if line 4 says "No data for your input selection" there was an
% error.
newlines = find(result == sprintf('\n'));
if length(newlines) < 71
    % Output the error with the HTML tags removed.
    error('iri:parseresult:modelWebError', regexprep(['The online ' ...
        'interface returned:\n' result], '<[^>]*>', ''));
elseif strfind( result( newlines(3)+1 : newlines(4) ), ...
        'No data for your input selection' )
    error('iri:parseresult:modelWebError', ['The online interface ' ...
         'returned: No data for your input selection.']);
end
data = sscanf(result(newlines(70)+1 : newlines(70 + sweeplen)-1), '%f');