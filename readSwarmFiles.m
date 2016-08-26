function swarmData = readSwarmFiles()

swarmFiles = dir('swarm/SW*');

lat = [];
lon = [];
lst = [];
alt = [];
rho = [];
timestamps = [];

targetCount = round(length(swarmFiles) / 10);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Reading SWARM, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );

parfor i = 1:length(swarmFiles)
    f = cdfread(['swarm/', swarmFiles(i).name]);
    N = size(f,1);
    tempLat = zeros(N,1);
    tempLon = zeros(N,1);
    tempLst = zeros(N,1);
    tempAlt = zeros(N,1);
    tempRho = zeros(N,1);
    tempTimestamps = zeros(N,1);
    for k = 1:N
        densThisRow = f{k,6};
        if densThisRow > 1E30
            continue;
        end
        tempTimestamps(k) = todatenum(f{k,1});
        tempAlt(k) = f{k,2}/1000;
        tempLat(k) = f{k,3};
        tempLon(k) = f{k,4};
        tempLst(k) = f{k,5};
        tempRho(k) = f{k,6};
    end
    lat = [lat; tempLat];
    lon = [lon; tempLon];
    lst = [lst; tempLst];
    alt = [alt; tempAlt];
    rho = [rho; tempRho];
    timestamps = [timestamps; tempTimestamps];
    
    if mod(i,10) == 0
        p.progress;
    end
end
p.stop;

[timestamps,order] = unique(timestamps);
lat = lat(order);
lon = lon(order);
lst = lst(order);
alt = alt(order);
rho = rho(order);

swarmData = struct('density', rho, 'latitude', lat, 'longitude', lon, 'solarTime', lst,...
                   'altitude', alt, 'timestamps', timestamps, 'densityError', zeros(size(rho)));

end