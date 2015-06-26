function [] = changeTIEGCMdate(datestring, inputFile)

% Check inputs
if ~(ischar(datestring) && ischar(inputFile))
    error('datestring and inputFile must all be strings')
end

filename = ['tiegcm_init_', datestring, '.nc'];
if exist(filename, 'file')
    delete(filename);
end
copyfile(inputFile, filename);

% Current datenum.
thisDatenum = datenum(datestring);
[~,~,~,hours,minutes,~] = datevec(thisDatenum);
if ~(hours == 0 && minutes == 0)
    error('Not allowed to specify hours, only date, e.g. 2010-12-31')
end

year = str2double(datestr(thisDatenum, 'yyyy'));
doy = floor(thisDatenum) - datenum(datestr(thisDatenum, 'yyyy'), 'yyyy') + 1;
[~,~,~,hours,minutes,~] = datevec(thisDatenum);
timestep = 30;

ncwrite(filename, 'year', year);
ncwrite(filename, 'timestep', timestep);
ncwrite(filename, 'ut', 0);
ncwrite(filename, 'time', 0);
modelTime = [doy-1; hours; minutes];
ncwrite(filename, 'mtime', modelTime);
iter = (modelTime(1) * 86400 + modelTime(2) * 3600 + modelTime(3) * 60) / timestep;
ncwrite(filename, 'iter', iter);
ncwrite(filename, 'day', modelTime(1));
ncwrite(filename, 'calendar_advance', 0);

end