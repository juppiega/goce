function [AE, AU, AL] = computeAE (filename)

fid = fopen('international_quiet_days.dat');
nlines = linecount(fid) - 1;
yr = zeros(nlines,1);
mo = zeros(nlines,1);
q = zeros(nlines,10);
tline = fgetl(fid);
tline = fgetl(fid);

k = 1;
while ischar(tline)
    yr(k) = str2double(tline(1:4));
    mo(k) = str2double(tline(6:7));
    q1 = str2double(tline(9:10));
    q2 = str2double(tline(11:12));
    q3 = str2double(tline(13:14));
    q4 = str2double(tline(15:16));
    q5 = str2double(tline(17:18));
    q6 = str2double(tline(20:21));
    q7 = str2double(tline(22:23));
    q8 = str2double(tline(24:25));
    q9 = str2double(tline(26:27));
    q10 = str2double(tline(28:29));
    q(k,:) = [q1,q2,q3,q4,q5,q6,q7,q8,q9,q10];
    tline = fgetl(fid);
    k = k + 1;
end
fclose(fid);

timestamps_quiet = datenum([yr,mo,ones(size(yr))]) + (q(:,1:5)-1);
timestamps_quiet = sort(timestamps_quiet(:));
timestamps_quiet_10 = datenum([yr,mo,ones(size(yr))]) + (q(:,1:10)-1);
timestamps_quiet_10 = sort(timestamps_quiet_10(:));

magfiles = dir(filename);
magFile = fopen(magfiles(1).name);
if magFile == -1
    error('supermag file open unsuccesful')
end
nlines = linecount(magFile) - 1;
N = 10000;
magdata = textscan(magFile, '%s %s %f32 %f32 %f32 %f32 %f32 %f32 %f32', N, 'Delimiter',',', 'HeaderLines',1);
X = nan(nlines,1,'single'); Y = X; station = X; mlat = X;
t_str = cell(nlines,1);
X(1:N) = single(magdata{7}); Y(1:N) = single(magdata{8}); mlat(1:N) = single(magdata{4}); t_str(1:N) = magdata{1}; station(1:N) = computeStation(magdata{2});
k = N;
while k + N < nlines 
    magdata = textscan(magFile, '%s %s %f32 %f32 %f32 %f32 %f32 %f32 %f32',N, 'Delimiter',',');
    X(k+1:k+N) = magdata{7}; Y(k+1:k+N) = magdata{8}; mlat(k+1:k+N) = single(magdata{4}); t_str(k+1:k+N) = magdata{1}; station(k+1:k+N) = computeStation(magdata{2});
    k = k+N;
end

r = nlines - k;
magdata = textscan(magFile, '%s %s %f32 %f32 %f32 %f32 %f32 %f32 %f32',r, 'Delimiter',',');
X(k+1:end) = magdata{7}; Y(k+1:end) = magdata{8}; mlat(k+1:end) = single(magdata{4}); t_str(k+1:end) = magdata{1}; station(k+1:end) = computeStation(magdata{2});
clear magdata

timestamps_all = datenum(t_str, 'yyyy-mm-dd HH:MM:SS'); clear t_str
timestamps = unique(timestamps_all);
stations = unique(station);

magarray = nan(length(timestamps), length(stations), 'single');

mlats = zeros(length(stations),1);
for k = 1:length(stations)
    ind = station == stations(k);
    H = single(sqrt(X(ind).^2 + Y(ind).^2));
    ind2 = ismember(timestamps, timestamps_all(ind));
    magarray(ind2, k) = H;
    mlats(k) = nanmedian(mlat(ind));
end


%stations = [];
%timestamps = [];
%for i = 1:length(magstruct)
%    names = fieldnames(magstruct(i));
%    stations = unique([stations; str2double(names)]);
%    t = [];
%    for k = 1:length(names)
%        t_this = getfield(magstruct(i), names{k}, 't');
%        t = unique([t; t_this]);
%    end
%    timestamps = unique([timestamps; t]);
%end
%
%magarray = nan(length(timestamps), length(stations));
%
%for i = 1:length(magstruct)
%    names = fieldnames(magstruct(i));
%    for k = 1:length(names)
%        t_this = getfield(magstruct(i), names{k}, 't');
%        H_this = getfield(magstruct(i), names{k}, 'H');
%        t_ind = ismember(timestamps, t_this);
%        H_ind = find(stations == str2double(names{k}));
%        magarray(t_ind, H_ind) = H_this;
%    end
%end

magarray = remove_baseline(magarray, timestamps, timestamps_quiet, timestamps_quiet_10);

%figure;
%plot(timestamps, magarray,'.')

N = length(timestamps);
AU = zeros(N,1);
AL = zeros(N,1);
AU_lat = AU;
AL_lat = AU;
num_stations = zeros(N,1);
for i = 1:length(timestamps)
    valid = ~isnan(magarray(i,:));
    if sum(valid) == 0
        warning(['Zero valid data for ', datestr(timestamps(i))])
    end
    [AU(i),k_AU] = nanmax(magarray(i,:));
    [AL(i),k_AL] = nanmin(magarray(i,:));
    num_stations(i) = sum(valid);
    AU_lat(i) = mlats(k_AU);
    AL_lat(i) = mlats(k_AL);
end

AE = AU - AL;

save([filename(1:end-4),'.mat'], 'timestamps','AE','AU','AL', 'num_stations', 'AU_lat', 'AL_lat')

%figure;
%plot(timestamps, [AE,AU,AL],'.')

end

function station = computeStation(station_string)

station = zeros(length(station_string),3);
for k = 1:length(station_string)
    station(k,:) = single(station_string{k});
end
station = single(10000*station(:,1) + 100*station(:,2) + station(:,3));

end

function magarray = remove_baseline(magarray, timestamps, timestamps_quiet, timestamps_quiet_10)

quiet_ind = false(size(timestamps));
quiet_ind_10 = false(size(timestamps));
for i = 1:length(timestamps_quiet)
    quiet_ind(timestamps_quiet(i) <= timestamps & timestamps < timestamps_quiet(i)+1) = true;
end
for i = 1:length(timestamps_quiet_10)
    quiet_ind_10(timestamps_quiet_10(i) <= timestamps & timestamps < timestamps_quiet_10(i)+1) = true;
end

[yr, mo, ~,~,~,~] = datevec(timestamps);
yrmo = yr*100+mo;
yrmo_un = unique(yrmo);
for i = 1:length(yrmo_un)
    t_ind = yrmo == yrmo_un(i);
    quiet_this = quiet_ind;
    quiet_this(~t_ind) = false;
    quiet_10 = quiet_ind_10;
    quiet_10(~t_ind) = false;
    quiet_mean = nan(size(magarray,2),1);
    for k = 1:size(magarray,2)
        H = magarray(t_ind, k);
        H_quiet = magarray(quiet_this, k);
        if sum(~isnan(H_quiet)) / length(H_quiet) < 0.6
            magarray(t_ind, k) = nan;
            continue;
        end
        H = H - mean(H_quiet(~isnan(H_quiet)));
        magarray(t_ind, k) = H;
        
        H_10 = magarray(quiet_10, k);
        quiet_mean(k) = mean(H_10(~isnan(H_10)).^2);
    end
    
    q_25 = quantile(quiet_mean(~isnan(quiet_mean)),0.25);
    q_75 = quantile(quiet_mean(~isnan(quiet_mean)),0.75);
    outliers = find(quiet_mean > q_75 + 1.0*(q_75-q_25));
    for k = 1:length(outliers)
        magarray(t_ind, outliers(k)) = nan;
    end
end

end
