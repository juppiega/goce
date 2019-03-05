function [AE, AU, AL] = computeAE (filename)

fid = fopen('international_quiet_days.dat');
nlines = linecount(fid) - 1;
yr = zeros(nlines,1);
mo = zeros(nlines,1);
q = zeros(nlines,5);
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
    q(k,:) = [q1,q2,q3,q4,q5];
    tline = fgetl(fid);
    k = k + 1;
end
fclose(fid);

timestamps_quiet = datenum([yr,mo,ones(size(yr))]) + (q-1);
timestamps_quiet = sort(timestamps_quiet(:));
timestamps = [];

magfiles = dir(filename);
magstruct = repmat(struct,length(magfiles),1);
for i = 1:length(magfiles)
    magFile = fopen(magfiles(i).name);
    if magFile == -1
        error('supermag file open unsuccesful')
    end
    magdata = textscan(magFile, '%s %s %f %f %f %f %f %f %f', 'Delimiter',',', 'HeaderLines',1);
    
    timestamps = datenum(magdata{1}, 'yyyy-mm-dd HH:MM:SS');
    station = zeros(length(timestamps),3);
    for k = 1:length(timestamps)
        station(k,:) = double(magdata{2}{k});
    end
    station = 1000000*station(:,1) + 1000*station(:,2) + station(:,3);
    stations = unique(station);
    for k = 1:length(stations)
        ind = station == stations(k);
        H = sqrt(magdata{7}(ind).^2 + magdata{8}(ind).^2);
        t = timestamps(ind);
        magstruct = setfield(magstruct, {i,1}, num2str(stations(k)), 'H', H);
        magstruct = setfield(magstruct, {i,1}, num2str(stations(k)), 't', t);
    end
    
end

stations = [];
timestamps = [];
for i = 1:length(magstruct)
    names = fieldnames(magstruct(i));
    stations = unique([stations; str2double(names)]);
    t = [];
    for k = 1:length(names)
        t_this = getfield(magstruct(i), names{k}, 't');
        t = unique([t; t_this]);
    end
    timestamps = unique([timestamps; t]);
end

magarray = nan(length(timestamps), length(stations));

for i = 1:length(magstruct)
    names = fieldnames(magstruct(i));
    for k = 1:length(names)
        t_this = getfield(magstruct(i), names{k}, 't');
        H_this = getfield(magstruct(i), names{k}, 'H');
        t_ind = ismember(timestamps, t_this);
        H_ind = find(stations == str2double(names{k}));
        magarray(t_ind, H_ind) = H_this;
    end
end

magarray = remove_baseline(magarray, timestamps, timestamps_quiet);

%figure;
%plot(timestamps, magarray,'.')

N = length(timestamps);
AU = zeros(N,1);
AL = zeros(N,1);
num_stations = zeros(N,1);
for i = 1:length(timestamps)
    valid = ~isnan(magarray(i,:));
    if sum(valid) == 0
        warning(['Zero valid data for ', datestr(timestamps(i))])
    end
    AU(i) = max(magarray(i,valid));
    AL(i) = min(magarray(i,valid));
    num_stations(i) = sum(valid);
end

AE = AU - AL;

save([filename(1:end-4),'.mat'], 'timestamps','AE','AU','AL', 'num_stations')

%figure;
%plot(timestamps, [AE,AU,AL],'.')

end

function magarray = remove_baseline(magarray, timestamps, timestamps_quiet)

quiet_ind = false(size(timestamps));
for i = 1:length(timestamps_quiet)
    quiet_ind(timestamps_quiet(i) <= timestamps & timestamps < timestamps_quiet(i)+1) = true;
end

[yr, mo, ~,~,~,~] = datevec(timestamps);
yrmo = yr*100+mo;
yrmo_un = unique(yrmo);
for i = 1:length(yrmo_un)
    t_ind = yrmo == yrmo_un(i);
    quiet_this = quiet_ind;
    quiet_this(~t_ind) = false;
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
        
        H_quiet = magarray(quiet_this, k);
        quiet_mean(k) = mean(abs(H_quiet(~isnan(H_quiet))));
    end
    
    q_25 = quantile(quiet_mean(~isnan(quiet_mean)),0.25);
    q_75 = quantile(quiet_mean(~isnan(quiet_mean)),0.75);
    outliers = find(quiet_mean > q_75 + 1.5*(q_75-q_25));
    for k = 1:length(outliers)
        magarray(t_ind, outliers(k)) = nan;
    end
end

end

function n = linecount(fid)

n = 0;
tline = fgetl(fid);
while ischar(tline)
  tline = fgetl(fid);
  n = n+1;
end
frewind(fid);

end
