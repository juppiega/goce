function out = iri2012(time, latitude, longitude, height, utc, coord, ...
    curldir, Rz12, IG12, f10_7_daily, f10_7_81day, tec_hmax, ne_top, ...
    fpeak, f2storm, bottom, f1_prob, auroral_boundary, foE_storm, ...
    Ne_Dreg, Te_top, ioncomp, NmF2_foF2, hmF2_M3000F2, NmE_foE, hmE)

% IRI2012 International Reference Ionosphere 2012 output parameters.
% 
% Usage: OUT = IRI2012(TIME, LATITUDE, LONGITUDE, HEIGHT, UTC, COORD,
%                  CURLDIR, Rz12, IG12, F10_7_DAILY, F10_7_81DAY, TEC_HMAX,
%                  Ne_TOP, FPEAK, F2STORM, BOTTOM, F1_PROB,
%                  AURORAL_BOUNDARY, foE_STORM, Ne_Dreg, Te_TOP, IONCOMP,
%                  NmF2_foF2, hmF2_M3000F2, NmE_foE, hmE)
% 
% Computes the International Reference Ionosphere (IRI) 2012, which is an
% internationally recognized model for various ionospheric properties. The
% position and time inputs can be scalars or arrays; when they are arrays,
% they should all have the same number of elements. The output is a matrix
% with 44 columns and the same number of rows as elements in the position
% and time inputs (possibly just 1). The ith row has the IRI output for the
% ith position/time input.
% 
% The function makes the computation by querying the online interface at
% http://omniweb.gsfc.nasa.gov/vitmo/iri2012_vitmo.html (hence internet
% access is required), which makes it pretty slow, especially when either
% more than one of the position and time inputs are made to vary or if the
% position and time inputs are not spaced linearly. If more than one input
% is varied, the function can be sped up by using a for loop and holding
% the smaller arrays constant, assuming the largest array is spaced
% linearly. This will result in fewer calls to the website since the
% website allows for a linear sweep in one variable. See the script
% iritest.m for an example of this.
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
%   determined by the input UTC. Valid range is from year 1958 to year 2016
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
%   -F10_7_DAILY: F10.7 radio flux, a daily index based on full disc solar
%   flux measurements at 2800 MHz (10.7 cm wavelength) first made at the
%   Algonquin Radio Observatory, near Ottawa, Canada (1947 until May 31,
%   1991) and then from the Dominion Radio Astrophysical Observatory, near
%   Penticton, British Columbia at local noon (1700 UT at Ottawa and 2000
%   UT at Penticton). The flux values are expressed in solar flux units 
%   (1 s.f.u. = 10-22 W*m-2*Hz-1 ). Indices are tabulated in two forms: the
%   "observed flux" (S), and the "adjusted flux" (Sa). The former are the
%   actual measured values, and are affected by the changing distance
%   between the Earth and Sun throughout the year, whereas the latter are
%   scaled to a standard distance of 1 AU. The "observed flux" values are
%   used in IRI. This index is typically calculated for you, but you can
%   enter your own and not rely on the internal index file. Valid range is
%   0 to 400.
%   
%   -F10_7_81DAY: The 81-day running mean of the daily F10.7 index
%   described above. 81 days are about 3 solar rotations and statistical
%   studies have found good correlation between this parameter or
%   PF10.7=(F10.7_daily + F10.7_81day)/2 and ionospheric parameters. This
%   index is typically calculated for you, but you can enter your own and
%   not rely on the internal index file. Valid range is 0 to 400.
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
%     3. ABT-2009: Altadill et al. [11] used a large volume of ionosonde
%     data to develop a much improved representation of latitudinal and
%     solar cycle variation of B0 and B1. B1 is a parameter describing the
%     bottomside profile shape. Selected by inputting one of 3, 'ABT-2009',
%     or 'ABT'.
%   The default option is ABT-2009.
% 
%   -F1_PROB: This parameter describes the occurrence probability of an F1
%   layer. There are three options here:
%     1. IRI-95: This option uses the ITU-recommended Ducharme et al. model
%     [12] that applies a simple cutoff solar zenith angle, so probability
%     is either 0 or 1. Selected by inputting one of 3, 'IRI-95', 'IRI95',
%     'IRI', or '95'.
%     2. Scotto-1997 no L: This option uses the model developed by Scotto
%     et al. [13] using only ionograms with clear F1 layer presence and
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
%   -AURORAL_BOUNDARY: Turn auroral boundary on (true or 'on') or off
%   (false or 'off') (optional, default is off). Leaving this off will
%   produce faster results.
%   
%   -foE_STORM: Turn the E peak auroral storm model on (true or 'on') or
%   off (false or 'off') (optional, default is on). The model was developed
%   by Mertens et al. [14] and Fernandez et al. [15] and describes the
%   average storm behavior in terms of the ratio foE_storm/foE_quite based
%   on the ap index history.
% 
%   -Ne_Dreg: D-region electron density. There are two options here:
%     1. FT-2001: Model based on Friedrich's compilation of reliable rocket
%     data [16]. Selected by inputting one of 2, 'FPT-2000', 'FPT', or
%     '2000'.
%     2. IRI-95: Model based on a much smaller selection of typical rocket
%     profiles [17]. Selected by inputting one of 1, 'IRI-95', 'IRI95',
%     'IRI', or '95'.
%   The default option is IRI-95.
% 
%   -Te_TOP: Topside electron temperature. There are two options here:
%     1. IRI-95 is using the global models at fixed altitudes developed by
%     Brace and Theis [18] based on their ISIS and AE-C electron
%     temperature measurements. The model is described in [19] and also in
%     the IRI-90 report that is available as PDF document from the
%     references section of the IRI homepage. Selected by inputting one of
%     2, 'IRI-95', 'IRI95', 'IRI', or '95'.
%     2. TBT-2012 model is newer and uses a large volume of satellite in
%     situ measurements and includes variations with solar activity [20].
%     Selected by inputting one of 1, 'TBT-2012', 'TBT', or '2012'.
%   The default option is TBT-2012.
% 
%   -IONCOMP: Ion composition. There are two options here:
%     1. DS95/DY85: Uses the Danilov and Smirnova model [21] based on their
%     compilation of rocket data in the region below the F-peak and the
%     Danilov and Yaichnikov model [22] based on their compilation of
%     Russian high-altitude rocket measurements. Selected by inputting one
%     of 2, 'DS95/DY85', or 'DS95'.
%     2. RBV10/TTS03: Based on an adjustment of the ion composition from
%     the FLIP-model photochemistry to the IRI electron density profile in
%     region below the F peak [23]. In the region above the F-peak the
%     model of Triskova et al. [24] based on satellite ion mass
%     spectrometer measurements. Selected by inputting one of 1,
%     'RBV10/TTS03', or 'RBV10'.
%   The default option is RBV10/TTS03.
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
%   -NmE_foE: A measured value to update the IRI profile to actual
%   conditions. Only valid when varying HEIGHT linearly. If this number is
%   between 20 and 10^8, it represents the E peak density (NmE) in cm^-3.
%   If instead it is between 0.1 and 14, it represents the E plasma
%   frequency in MHz.
% 
%   -hmE: Measured E peak height in km to update the IRI profile to actual
%   conditions. Valid range is between 70 km and 200 km.
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
%     34. F1 probability.
%     35. CGM latitude of auroral oval boundary.
%     36. Spread F probability.
%     37. Ratio of foE storm to foE quiet.
%     38. 12-month running mean of sunspot number Rz12 used by the model.
%     39. Ionospheric index IG12 used by the model.
%     40. Daily solar radio flux F107D used by the model.
%     41. 81 day solar radio flux F107_81D used by the model.
%     42. 3 hour ap index used by the model.
%     43. Daily ap index used by the model.
%     44. 3 hour Kp index used by the model.
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
%   [11] D. Altadilla, J. M. Tortaa, and E. Blanch, "Proposal of new models
%   of the bottom-side B0 and B1 parameters for IRI,", _Adv._Space_Res._,
%   vol. 43, no. 11, 2009, doi:10.1016/j.asr.2008.08.014.
%   [12] E. D. DuCharme, L. E. Petrie, and R. Eyfrig, "A method for
%   predicting the F1 layer critical frequency based on the Zurich smoothed
%   sunspot number," _Radio_Sci._, vol. 8, no. 10, pp. 837-839, 1973,
%   doi:10.1029/RS008i010p00837.
%   [13] C. Scotto, S. M. Radicella, and B. Zolesi, "An improved
%   probability function to predict the F1 layer occurrence and L
%   condition," _Radio_Sci._, vol. 33, no. 6, pp. 1763-1765, 1998,
%   doi:10.1029/98RS02637.
%   [14] C. J. Mertens, X. Xu, D. Bilitza, M. G. Mlynczak, and J. M. Russel
%   III, "Empirical STORM-E model: I. Theoretical and observational basis,"
%   _Adv._Space_Res._, vol 51, no. 4, pp. 554-574,
%   doi:10.1016/j.asr.2012.09.009.
%   and C. J. Mertens, X. Xu, D. Bilitza, M. G. Mlynczak, and J. M. Russel
%   III, "Empirical STORM-E model: II.  Geomagnetic corrections to
%   nighttime ionospheric E-region electron densities," _Adv._Space_Res._,
%   vol 51, no. 4, pp. 575-98, doi:10.1016/j.asr.2012.09.014.
%   [15] J. R. Fernandez, C. J. Mertens, D. Bilitza, X. Xu, J. M. Russel
%   III, M. G. Mlynczak, "Feasibility of developing an ionospheric E-region
%   electron density storm model using TIMED/SABER measurements," vol 46,
%   no. 8, pp. 1070-1077, doi:10.1016/j.asr.2010.06.008.
%   [16] M. Friedrich and K. M. Torkar, "FIRI: A semiempirical model of the
%   lower ionosphere," _J._Geophys._Res._, vol 106, no. A10, 2001,
%   doi:10.1029/2001JA900070.
%   [17] E. Mechtly and D. Bilitza, "Models of D-region electron
%   concentration," Rep. IPW-WB1, Inst. fur phys. Weltraumforsch.,
%   Freiburg, Germany, 1974.
%   [18] L. H. Brace and R. F. Theis, "Global empirical models of
%   ionospheric electron temperature in the upper F region and plasmasphere
%   based on in situ measurements from the Atmospheric Explorer-C, ISIS-1
%   and ISIS-2 satellites," _J._Atmos._Terr._Phys._, vol. 43, pp.
%   1317–1343, 1981.
%   [19] D. Bilitza, L. H. Brace, R. F. Theis, "Modelling of ionospheric
%   temperature profiles," _Adv._Space_Res._, vol. 5, no. 7, pp. 53-58,
%   1985, doi:10.1016/0273-1177(85)90356-4.
%   [20] V. Truhlík, D. Bilitza, and L. T?ísková, "A new global empirical
%   model of the electron temperature with the inclusion of the solar
%   activity variations for IRI," _Earth_Planets_Space_, vol. 64, no. 6,
%   pp. 531-543, 2012 doi:10.5047/eps.2011.10.016.
%   [21] A. D. Danilov, and N. V. Smirnova, "Improving the 75 to 300 km ion
%   composition model of the IRI," _Adv._Space_Res._, vol. 15, no. 2, pp.
%   171-177, 1995, doi:10.1016/S0273-1177(99)80044-1.
%   [22] A. D. Danilov, and A. P. Yaichnikov, "A new model of the ion
%   composition model at 75 to 1000 km for IRI," _Adv._Space_Res._, 
%   vol. 5, no. 7, pp. 75–79, 1985, doi:10.1016/0273-1177(85)90360-6.
%   [23] P. G. Richards, D. Bilitza, and D. Voglozin, "Ion density
%   calculator (IDC): A new efficient model of ionospheric ion densities,"
%   _Radio_Sci._, vol. 45, p. RS5007, 2010, doi:10.1029/2009RS004332.
%   [24] V. Truhlík, L. T?ísková, J. Šmilauer, "An empirical model of
%   ion composition in the outer ionosphere," _Adv._Space_Res._, vol. 31,
%   no. 3, pp. 653-663, doi:10.1016/S0273-1177(03)00040-1.
% 
% See also: IRI2007, IRITEST, MSIS, IGRF, DATENUM.

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
    f10_7_daily = [];
end
if nargin < 11
    f10_7_81day = [];
end
if nargin < 12
    tec_hmax = [];
end
if nargin < 13 || isempty(ne_top)
    ne_top = 'NeQuick';
end
if nargin < 14 || isempty(fpeak)
    fpeak = 'URSI';
end
if nargin < 15 || isempty(f2storm)
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
if nargin < 16 || isempty(bottom)
    bottom = 'ABT-2009';
end
if nargin < 17 || isempty(f1_prob)
    f1_prob = 'Scotto-1997 no L';
end
if nargin < 18 || isempty(auroral_boundary)
    auroral_boundary = false;
elseif ischar(auroral_boundary)
    if strcmpi(auroral_boundary, 'on')
        auroral_boundary = true;
    elseif strcmpi(auroral_boundary, 'off')
        auroral_boundary = false;
    else
        error('iri:invalid_auroral_boundary', ...
            'Unrecognized input ''%s'' for AURORAL_BOUNDARY.', ...
            auroral_boundary);
    end
end
if nargin < 19 || isempty(foE_storm)
    foE_storm = true;
elseif ischar(foE_storm)
    if strcmpi(foE_storm, 'on')
        foE_storm = true;
    elseif strcmpi(foE_storm, 'off')
        foE_storm = false;
    else
        error('iri:invalid_foE_storm', ...
            'Unrecognized input ''%s'' for foE_STORM.', foE_storm);
    end
end
if nargin < 20 || isempty(Ne_Dreg)
    Ne_Dreg = 'IRI-95';
end
if nargin < 21 || isempty(Te_top)
    Te_top = 'TBT-2012';
end
if nargin < 22 || isempty(ioncomp)
    ioncomp = 'RBV10/TTS03';
end
if nargin < 23 || isempty(NmF2_foF2)
    NmF2_foF2 = 0;
end
if nargin < 24 || isempty(hmF2_M3000F2)
    hmF2_M3000F2 = 0;
end
if nargin < 25 || isempty(NmE_foE)
    NmE_foE = 0;
end
if nargin < 26 || isempty(hmE)
    hmE = 0;
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
if isempty(f10_7_daily) || (f10_7_daily >= 0 && f10_7_daily <= 400)
    f10_7_daily = sprintf('&radio_f=%#g', f10_7_daily);
else
    error('iri:invalidf10_7_daily', ...
        'Input f10_7_daily is %g must be between 0 and 400.', f10_7_daily);
end
if isempty(f10_7_81day) || (f10_7_81day >= 0 && f10_7_81day <= 400)
    f10_7_81day = sprintf('&radio_f81=%#g', f10_7_81day);
else
    error('iri:invalidf10_7_81day', ...
        'Input f10_7_81day is %g must be between 0 and 400.', f10_7_81day);
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
    case {3, 'abt-2009', 'abt'}
        bottom = '2.';
    otherwise
        if ~ischar(bottom)
            bottom = num2str(bottom);
        else
            bottom = ['' bottom ''];
        end
        error('iri:invalidBottom', ['Unrecognized command ' bottom ...
            ' for input BOTTOM. Valid options are B0 Table (1, ' ...
            '''B0 Table'', ''B0Table'', ''B0'', or ''B''), Gulyaeva ' ...
            '(2, ''Gulyaeva'', ''Guly'', or ''G''), and ABT-2009 ' ...
            '(3, ''ABT-2009'', or ''ABT'').']);
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
if auroral_boundary
    auroral_boundary = '0.';
else
    auroral_boundary = '1.';
end
if foE_storm
    foE_storm = '0.';
else
    foE_storm = '1.';
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
    case {1, 'tbt-2012', 'tbt', '2012'}
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
            ' for input TE_TOP. Valid options are TBT-2012 (1, ' ...
            '''TBT-2012'', ''TBT'', or ''2012'') and IRI-95 (2, ' ...
            '''IRI-95'', ''IRI95'', ''IRI'', or ''95'').']);
end
switch lower(ioncomp)
    case {1, 'rbv10/tts03', 'rbv10'}
        ioncomp = '0.';
    case {2, 'ds95/dy85', 'ds95'}
        ioncomp = '1.';
    otherwise
        if ~ischar(ioncomp)
            ioncomp = num2str(ioncomp);
        else
            ioncomp = ['' ioncomp ''];
        end
        error('iri:invalidIoncomp', ['Unrecognized command ' ioncomp ...
            ' for input IONCOMP. Valid options are RBV10/TTS03 (1, ' ...
            '''RVB10/TTS03'', ''RBV10'') and DS95/DY85 (2, ' ...
            '''DS95/DY85'' or ''DS95'').']);
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
if ~(NmE_foE == 0 || (NmE_foE >= 20 && NmE_foE <= 1e8) || ...
        (NmE_foE >= 0.1 && NmE_foE <= 14))
    error('iri:invalidNmE_foE', ['Input NmE_foE is %g but must be ' ...
        'one of 0, between 20 and 10^8 (cm^-3), or between 0.1 and 14 ' ...
        '(MHz).'], ...
        NmE_foE);
end
if ~(hmE == 0 || (hmE >= 70 && hmE <= 200))
    error('iri:invalidhmE', ['Input hmE is %g but must be either 0 or ' ...
        'between 70 and 200 (km).'], hmE);
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
endcmd = [Rz12 IG12 f10_7_daily f10_7_81day tec_hmax ...
    '&ne_top=' ne_top ...
    '&imap=' fpeak ...
    '&ffof2=' f2storm ...
    '&ib0=' bottom ...
    '&probab=' f1_prob ...
    '&fauroralb=' auroral_boundary ...
    '&ffoE=' foE_storm ...
    '&dreg=' Ne_Dreg ...
    '&tset=' Te_top ...
    '&icomp=' ioncomp ...
    sprintf('&nmf2=%#g', NmF2_foF2) ...
    sprintf('&f2=%#g', hmF2_M3000F2) ...
    sprintf('&user_nme=%#g', NmE_foE) ...
    sprintf('&user_hme=%#g', hmE) ...
    sprintf('&vars=%i', 17:60) ...11:50) ...
    '" http://omniweb.gsfc.nasa.gov/cgi/vitmo/vitmo_model.cgi > "' ...
    ...'" http://ccmc.gsfc.nasa.gov/cgi-bin/modelweb/models/vitmo_model.cgi > "' ...
    fname '"'];

% Get year, month, day, hour for IRI.
[year, month, day, hour, minute, second] = datevec(time);
hour = hour + minute/60 + second/3600;
dayyear = dayofyear(year, month, day);

% Get the beginning of the string to be input into the system function.
if isunix || ismac % Curl is built-in to operating system.
    initialcmd = 'curl -d "model=iri_2012';
else % Use the input curldir (possibly the default for that input).
    initialcmd = ['"' fullfile(curldir, 'curl"') ' -d "model=iri_2012'];
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
                'It returned status:%i, %s.'], status, result);
        end
    case {1 2 3 4 5 7 8}
        % The online interface will only output up to a sweep length of
        % 1000, so split the sweep into increments of 1000.
        nsweeps = ceil(numel(sweep) / 1e3);
        sweepstart = sweep(1 : 1e3 : end);
        sweepstop = sweep(1e3 : 1e3 : end);
        if numel(sweepstart) - 1 == numel(sweepstop)
            sweepstop(end + 1) = sweep(end);
        end
        sweepstep = mode(diff(sweep));
        sweeplen = round((sweepstop - sweepstart) ./ sweepstep + 1);
        prevsweeplen = [0 sweeplen(1:end-1)];
        sweepstop = min([sweep(end) + sweepstep/10, sweepmax]);
        out = zeros(44*sum(sweeplen), 1);
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
                out((1:44*sweeplen(index)) + 44*prevsweeplen(index)*...
                    (index-1)) = parseresult(fname, sweeplen(index));
            else
                error('iri:curlError', ['Curl command did not work. ' ...
                    'It returned status:%i, %s.'], status, result);
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
        out = zeros(44*maxnum, 1);
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
                out((1:44) + 44*(index-1)) = parseresult(fname, 1);
            else
                error('iri:curlError', ['Curl command did not work. ' ...
                    'It returned status:%i, %s.'], status, result);
            end
        end
end % End case

% Right now, out has all the data in one long vector, but we want it to be
% an array, so reshape.
out = reshape(out, 44, []).';

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

% Data starts at line 78. If there isn't a line 78, there was an error.
% Also, if line 4 says "No data for your input selection" there was an
% error.
newlines = find(result == sprintf('\n'));
if length(newlines) < 78
    % Output the error with the HTML tags removed.
    error('iri:parseresult:modelWebError', regexprep(['The online ' ...
        'interface returned:\n' result], '<[^>]*>', ''));
elseif strfind( result( newlines(3)+1 : newlines(4)-1 ), ...
        'No data for your input selection' )
    error('iri:parseresult:modelWebError', ['The online interface ' ...
         'returned: No data for your input selection.']);
end
data = sscanf(result(newlines(77)+1 : newlines(77 + sweeplen)-1), '%f');