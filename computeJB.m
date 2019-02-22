function rho = computeJB(S)

poolobj = gcp('nocreate'); % If no pool, create a new one.
if isempty(poolobj)
    parpool();
end

[yr,~,~,h,m,s] = datevec(S.timestamps);
% Annual parameter.
if ~isfield(S, 'doy')
    yearVec = [yr, repmat([1,1,0,0,0], length(yr), 1)];
    S.doy = S.timestamps - datenum(yearVec) + 1;
end

N = length(S.latitude);
rho = zeros(N,1);

targetCount = round(N / 10000);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running JB2008, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );
                
doy = S.doy;
alt = S.altitude;
lat = S.latitude;
lon = S.longitude;
parfor i = 1:N
    [~,~,~,~,~,~,~,~,~,~,~,~,rho(i),~,~,~,~,~,~,~,~,~,~] = jb2k8(yr(i),doy(i),h(i),m(i),s(i),alt(i),lat(i),lon(i));
    if mod(i, 10000) == 0
        p.progress;
    end
end
p.stop;



end