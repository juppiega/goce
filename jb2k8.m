function [outputyear,outputdoy,outputhr,outputmin,outputsec,outputht,outputlat,outputlon,outputsun1,outputsun2,outputtemp1,outputtemp2,outputdensity,outputD1950,outputF10,outputF10B,outputS10,outputS10B,outputXM10,outputXM10B,outputY10,outputY10B,outputDSTDTC] = jb2k8(years,doy,hour,min,sec,ht,lat,lon)

% JB2008 Calling Script.
% jb2k8 is used to call the JB2008 model in MATLAB. It can also be embedded
% in a SIMULINK block if desired.
%
% Replace any outputs with a ~ to skip them.
% [~,~,~,~,~,~,~,~,~,~,~,~,outputdensity,~,~,~,~,~,~,~,~,~,~] = jb2k8(years,doy,hour,min,sec,ht,lat,lon)
%
%
%



%% Check for latest version and all required files.
% The sum of the "exist" results should add to 8.

FileSum = exist('SOLFSMY.TXT') + exist('SOLRESAP.TXT') + ...
 exist('DSTFILE.TXT') +  exist('DTCFILE.TXT');

if  verLessThan('matlab', '8.4') == 1 & FileSum == 8
    disp(['MATLAB 2014b or later is required to automatically update files.'])
    disp(['Files may be out of date!'])
end


%% Update required text files if necessary

UpdateFile('SOLFSMY.TXT');
UpdateFile('DTCFILE.TXT');
UpdateFile('DSTFILE.TXT');
UpdateFile('SOLRESAP.TXT');

%% Check the the input date is not later than the latest date of the input files.

doytoday = datenum(date) - datenum(year(date),1,1) + 1;

if years >= year(date) & doy >= doytoday
    disp(['Choose another date.'])
return
end

%% Call the JB2008 function

formatSpec = './JB2008%5d%5d%5d%5d%7.3f%9.2f%8.2f%8.2f';
inputdata = sprintf(formatSpec,years,floor(doy),hour,min,sec,ht,lat,lon);

[status,cmdout] = system(inputdata);
if status
    disp(cmdout)
    error('JB2008 command line fail')
end

G = textscan(cmdout,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');


%% Outputs
outputs1 = G{1,1};
outputs2 = G{1,2};
if length(outputs1) < 21
    disp(inputdata)
    disp(cmdout)
    disp(outputs1)
    error('JB2008 output too short')
end

% Time
    outputyear = outputs1(1);
    outputdoy = outputs1(2);
    outputhr = outputs1(3);
    outputmin = outputs1(4);
    outputsec = outputs1(5);
% Location
    outputht = outputs1(6);
    outputlat = outputs1(7);
    outputlon = outputs1(8);
% Misc
    outputsun1 = outputs1(9);
    outputsun2 = outputs2(9);
    outputtemp1 = outputs1(10);
    outputtemp2 = outputs2(10);

% Density
    outputdensity  = outputs1(11);

% Indices
    outputD1950 = outputs1(12);
    outputF10 = outputs1(13);
    outputF10B = outputs1(14);
    outputS10 = outputs1(15);
    outputS10B = outputs1(16);
    outputXM10 = outputs1(17);
    outputXM10B = outputs1(18);
    outputY10 = outputs1(19);
    outputY10B = outputs1(20);
    outputDSTDTC = outputs1(21);
    
end

function [y] = year(date)

vec = datevec(date);
y = vec(1);

end

