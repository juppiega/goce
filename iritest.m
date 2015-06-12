% IRITEST Plot electron density around the world to test IRI functions.

%% Clear variables and close all figures.
tic;
clear;
close all;

%% Test values.
time = datenum([2012 7 17 12 0 0]); % Must be a scalar in this script.
latitude = -90:10:90;    % Degrees. Can be a vector or scalar.
longitude = -180:10:180; % Degrees. Can be a vector or scalar.
altitude = 100;         % km. Must be a scalar in this script.
fun2test = @iri2012; % Either @iri2007 or @iri2012.
useparfor = true; % True: use parfor (requires Parallel Computing Toolbox).
                  % False: use a regular for loop.

% %% Compute B-field (needs function IGRF from MATLAB file exchange).
% [LAT, LON] = meshgrid(latitude, longitude); s = size(LAT);
% [Bx, By, Bz] = igrf(time, LAT, LON, altitude, 'geod');
% Bx = reshape(Bx, s); By = reshape(By, s); Bz = reshape(Bz, s);
% decl = atan(By./Bx)*180/pi;
% B = hypot(Bx, hypot(By, Bz));

%% Compute electron density for all values above.
Ne = zeros(numel(latitude), numel(longitude));

% Use a parfor loop or regular for loop?
if useparfor && license('test', 'Distrib_Computing_Toolbox')
    
    % Check for the parfor_progress function.
    use_parfor_progress = exist('parfor_progress.m', 'file');
    
    % Make the function call. Index into the vector (latitude or longitude)
    % with fewer elements.
    if numel(latitude) >= numel(longitude)% || altitude ~= 100
        latitude = latitude(:); N = numel(longitude);
        if use_parfor_progress
            temp = parfor_progress(N); clear temp; %#ok
        end
        parfor index = 1:N
            out = fun2test(time, latitude, longitude(index), altitude);
            Ne(:, index) = out(:, 1);
            if use_parfor_progress
                fprintf('%4.4g%% Complete\n', parfor_progress);
            end
        end
    else
        longitude = longitude(:); N = numel(latitude);
        if use_parfor_progress
            temp = parfor_progress(N); clear temp; %#ok
        end
        parfor index = 1:N
            out = fun2test(time, latitude(index), longitude, altitude);
            Ne(index, :) = out(:, 1).';
            if use_parfor_progress
                fprintf('%4.4g%% Complete\n', parfor_progress);
            end
        end
    end
    
    % Close the parfor_progress.txt file.
    if use_parfor_progress
        temp = parfor_progress(0); clear temp; %#ok
    end
    
else % Use regular for loop.
    
    % Index into the vector (latitude or longitude) with fewer elements.
    if numel(latitude) >= numel(longitude)% || altitude ~= 100
        latitude = latitude(:); N = numel(longitude);
        for index = 1:N
            out = fun2test(time, latitude, longitude(index), altitude);
            Ne(:, index) = out(:, 1);
            fprintf('%4.4g%% Complete\n', index/N*100);
        end
    else
        longitude = longitude(:); N = numel(latitude);
        for index = 1:N
            out = fun2test(time, latitude(index), longitude, altitude);
            Ne(index, :) = out(:, 1).';
            fprintf('%4.4g%% Complete\n', index/N*100);
        end
    end
    
end

%% Plot data.
figure;
if ~license('test', 'MAP_Toolbox')
    hold on;
    % surf(LON.', LAT.', B.'); %decl
    surf(longitude, latitude, Ne.');
    shading flat;
    clim = get(gca, 'CLim'); zlim = get(gca, 'ZLim');
    load('topo.mat', 'topo'); topo = [topo(:, 181:360), topo(:, 1:180)];
    [C, h] = contour3(-180:179, -90:+89, topo + zlim(2), [0 0] + zlim(2));
    set(h, 'EdgeColor', 0.25*[1 1 1]);
    set(gca, 'XLim', [-180 180], 'XTick', [], ...
        'YLim', [-90 90], 'YTick', [], 'CLim', clim);
else
    axesm miller; axis fill;
    hold on;
    % surfm(LAT, LON, B); %decl
    surfm(latitude, longitude, Ne);
    load coast;
    plotm(lat, long, 'Color', 0.25*[1 1 1]);
end
hc = colorbar;
title(hc, '\itN_e\rm in m^{-3}');
% title(['Magnetic Field Declination in degrees at \ith\rm = ' ...
%     sprintf('%g km', altitude) ' at ' datestr(time) ' UTC']);
title(['Electron Density (\itN_e\rm) at \ith\rm = ' ...
    sprintf('%g km', altitude) ' at ' datestr(time) ' UTC']);
print -dpng -r100 iri.png

%%
toc;