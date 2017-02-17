function computeT0()

z0 = 130;

files = [dir('EISCAT_*'),dir('MILLSTONE_HILL*'),dir('ARECIBO*')];
for i = 1:length(files)
    data = load(files(i).name);
    if size(data,2) ~= 13
        warning(['File: ', files(i).name,' did not contain any data!'])
        continue;
    end
    year = data(:,1);
    month = data(:,2);
    day = data(:,3);
    hour = data(:,4);
    minutes = data(:,5);
    seconds = data(:,6);
    altitude = data(:,7);
    vertRes = data(:,8);
    fundPulsLen = data(:,9);
    Ti = data(:,10);
    Ti_err = data(:,11);
    latitude = data(:,12);
    longitude = data(:,13);
    
    timestamps = datevec([year, month, day, hour, minutes, seconds]);
    solarTime = hour + longitude/15 + minutes/60 + seconds/3600;
    solarTime(solarTime >= 24) = solarTime(solarTime >= 24) - 24;
    solarTime(solarTime < 0) = solarTime(solarTime < 0) + 24;
    
    justBelowZ0 = find((diff(double(altitude > z0)) ~= 0));
    justAboveZ0 = justBelowZ0 + 1;
    if justAboveZ0(end) > length(altitude)
        justAboveZ0(end) = [];
        justBelowZ0(end) = [];
    end
    
    TiJustBelow = Ti(justBelowZ0);
    TiJustAbove = Ti(justAboveZ0);
    TiErrJustBelow = Ti_err(justBelowZ0);
    TiErrJustAbove = Ti_err(justAboveZ0);
    
    rmFromInterp = (isnan(TiJustBelow) | isnan(TiJustAbove) | ...
                    isnan(TiErrJustBelow) | isnan(TiErrJustAbove));
    justBelowZ0(rmFromInterp) = [];
    justAboveZ0(rmFromInterp) = [];
    
    TiJustBelow = Ti(justBelowZ0);
    TiJustAbove = Ti(justAboveZ0);
    TiErrJustBelow = Ti_err(justBelowZ0);
    TiErrJustAbove = Ti_err(justAboveZ0);
    altJustBelow = altitude(justBelowZ0);
    altJustAbove = altitude(justAboveZ0);
    
    solarTime = solarTime(justBelowZ0);
    timestamps = timestamps(justBelowZ0);
    longitude = longitude(justBelowZ0);
    latitude = latitude(justBelowZ0);
    
    
    N = size(justBelowZ0);
    Ti = zeros(N);
    Ti_err = zeros(N);
    altitude = zeros(N) + z0;
    
    m = (TiJustAbove - TiJustBelow) ./ (altJustAbove - altJustBelow);
    Ti = m.*(z0 - altJustBelow) + TiJustBelow;
    Ti_err = sqrt(TiErrJustAbove.^2 + TiErrJustBelow.^2);
    vertRes = max(vertRes(justBelowZ0), vertRes(justAboveZ0));
    fundPulsLen = max(fundPulsLen(justBelowZ0), fundPulsLen(justAboveZ0));
    
    if ~isempty(strfind(files(i).name, 'EISCAT'))
        index = zeros(N) + 1;
    elseif ~isempty(strfind(files(i).name, 'MILLSTONE'))
        index = zeros(N) + 2;
    elseif ~isempty(strfind(files(i).name, 'ARECIBO'))
        index = zeros(N) + 3;
    end
    
    T0(i) = struct('data', Ti, 'Ti', Ti, 'Ti_err', Ti_err, 'Tn_err', Ti_err,...
        'solarTime', solarTime, 'timestamps', timestamps, 'longitude', longitude,...
        'latitude', latitude, 'altitude', altitude, 'pulseLen', fundPulsLen, ...
        'fundPulsLen', fundPulsLen, 'vertRes', vertRes, 'index', index);
    
end

save('T0.mat', 'T0')

end