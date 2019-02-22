function rho = computeJB(S)

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
                

for i = 1:N
    [~,~,~,~,~,~,~,~,~,~,~,~,rho(i),~,~,~,~,~,~,~,~,~,~] = jb2k8(yr(i),S.doy(i),h(i),m(i),s(i),S.altitude(i),S.latitude(i),S.longitude(i));
    if mod(i, 10000) == 0
        p.progress;
    end
    i
end
p.stop;



end