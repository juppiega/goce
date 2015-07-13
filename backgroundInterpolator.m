function [] = backgroundInterpolator(filenames)

tiegcmOldFiles = {};
load('goceVariables.mat', 'latitude', 'longitude', 'altitude', 'timestampsDensityDatenum');

if strcmpi(pause('query'), 'off')
    pause on;
end

while 1
    tiegcmCurrentFiles = dir(filenames);
    newFiles = findNewFiles(tiegcmCurrentFiles, tiegcmOldFiles);
    if any(newFiles)
        if exist('tiegcmDens.mat', 'file')
            load tiegcmDens
        else
            tiegcmGoceDatenums = []; tiegcmGoce270km = []; tiegcmGoceInterp = [];
            save('tiegcmDens.mat', 'tiegcmGoceDatenums', '-v7')
            save('tiegcmDens.mat', 'tiegcmGoce270km', '-append')
            save('tiegcmDens.mat', 'tiegcmGoceInterp', '-append')
            load tiegcmDens
        end

        newFileNames = tiegcmCurrentFiles(newFiles);
        lat = ncread(newFileNames(1).name, 'lat');
        lon = ncread(newFileNames(1).name, 'lon');
        for i = 1:length(newFileNames)
            thisFileDatenums = giveTiegcmDatenums(newFileNames(i).name);
            if interpolateThisFile(thisFileDatenums, tiegcmGoceDatenums, timestampsDensityDatenum)
                
                ind = (timestampsDensityDatenum >= thisFileDatenums(1) & timestampsDensityDatenum <= thisFileDatenums(end));
                goceDatenums = timestampsDensityDatenum(ind);
                goceLon = longitude(ind);
                goceLat = latitude(ind);
                goceAlt = altitude(ind);
                tiegcmAlt = double(ncread(newFileNames(i).name, 'ZG')) / 100;
                tiegcmDens = double(ncread(newFileNames(i).name, 'DEN'));
                tiegcmGoceInterpThisFile = nan(size(goceDatenums));
                tiegcmGoce270kmThisFile = nan(size(goceDatenums));
                warning('off', 'MATLAB:qhullmx:InternalWarning');
                fprintf('Interpolating new file: %s\n', newFileNames(i).name)
                %showInterval = round(length(goceDatenums) / 5);

                for j = 1:length(goceDatenums)
                    tiegcmGoceInterpThisFile(j) = interpSatellite(lon, lat, tiegcmAlt, thisFileDatenums, tiegcmDens,...
                                                            goceLon(j), goceLat(j), goceAlt(j), goceDatenums(j), 1);
                    tiegcmGoce270kmThisFile(j) = interpSatellite(lon, lat, tiegcmAlt, thisFileDatenums, tiegcmDens,...
                                                            goceLon(j), goceLat(j), 270E3, goceDatenums(j), 1);
                    if j == 1 || mod(j, 1000) == 0 || j == length(goceDatenums)
                        fprintf('%d %s %s\n', j, ' / ', num2str(length(goceDatenums)))
                    end
                end
                
                notNan = ~isnan(tiegcmGoceInterpThisFile);
                if length(find(notNan)) > 1
                    tgNoNans = tiegcmGoceInterpThisFile(notNan);
                    tNoNans = goceDatenums(notNan);
                    tiegcmGoceInterpThisFile = interp1(tNoNans, tgNoNans, goceDatenums, 'nearest', 'extrap');
                elseif length(find(notNan)) == 1 && length(goceDatenums) > 1
                    tgNoNan = tiegcmGoceInterpThisFile(notNan);
                    tiegcmGoceInterpThisFile(~notNan) = ones(length(~notNan),1) * tgNoNan;
                elseif length(find(notNan)) == 0
                    [~,nearInd] = min(abs(tiegcmGoceDatenums - goceDatenums));
                    tiegcmGoceInterpThisFile(~notNan) = ones(length(~notNan),1) * tiegcmGoceInterp(nearInd);
                elseif length(goceDatenums) == 0
                    continue;
                end
                
                notNan = ~isnan(tiegcmGoce270kmThisFile);
                if length(find(notNan)) > 1
                    tgNoNans = tiegcmGoce270kmThisFile(notNan);
                    tNoNans = goceDatenums(notNan);
                    tiegcmGoce270kmThisFile = interp1(tNoNans, tgNoNans, goceDatenums, 'nearest', 'extrap');
                elseif length(find(notNan)) == 1 && length(goceDatenums) > 1
                    tgNoNan = tiegcmGoce270kmThisFile(notNan);
                    tiegcmGoce270kmThisFile(~notNan) = ones(length(~notNan),1) * tgNoNan;
                elseif length(find(notNan)) == 0
                    [~,nearInd] = min(abs(tiegcmGoceDatenums - goceDatenums));
                    tiegcmGoce270kmThisFile(~notNan) = ones(length(~notNan),1) * tiegcmGoce270km(nearInd);
                elseif length(goceDatenums) == 0
                    continue;
                end
                
                tiegcmGoceDatenums = [tiegcmGoceDatenums; goceDatenums];
                tiegcmGoce270km = [tiegcmGoce270km; tiegcmGoce270kmThisFile];
                tiegcmGoceInterp = [tiegcmGoceInterp; tiegcmGoceInterpThisFile];

                [tiegcmGoceDatenums, uniqueInd] = unique(tiegcmGoceDatenums);
                tiegcmGoce270km = tiegcmGoce270km(uniqueInd);
                tiegcmGoceInterp = tiegcmGoceInterp(uniqueInd);
            end
        end
        save('tiegcmDens.mat','tiegcmGoceDatenums','-append')
        save('tiegcmDens.mat','tiegcmGoce270km','-append')
        save('tiegcmDens.mat','tiegcmGoceInterp','-append')
        tiegcmOldFiles = tiegcmCurrentFiles;
        copyfile('tiegcmDens.mat','backup_tiegcmDens.mat');
    end
    pause(5);
end

end

function interpolateThisFile = interpolateThisFile(thisFileDatenums, tiegcmDatenums, goceDatenums)

if isempty(tiegcmDatenums)
    interpolateThisFile = true;
    return
end

ind = (goceDatenums >= tiegcmDatenums(1) & goceDatenums <= tiegcmDatenums(end));
goceTimestamps = round((goceDatenums(ind) - tiegcmDatenums(1)) * 86400);
tiegcmTimestamps = round((tiegcmDatenums - tiegcmDatenums(1)) * 86400);
fileTimestamps = round((thisFileDatenums - tiegcmDatenums(1)) * 86400);

goceNotInterpolated = goceTimestamps(~ismember(goceTimestamps, tiegcmTimestamps));
if any(goceNotInterpolated >= fileTimestamps(1) & goceNotInterpolated <= fileTimestamps(end))
    interpolateThisFile = true;
end

end

function [newFiles] = findNewFiles(tiegcmCurrentFiles, tiegcmOldFiles)

oldFileNames = cell(length(tiegcmOldFiles), 1);
for i = 1:length(tiegcmOldFiles)
    oldFileNames{i} = tiegcmOldFiles(i).name;
end

newFiles = false(length(tiegcmCurrentFiles),1);
for i = 1:length(tiegcmCurrentFiles)
    oldsContainThisFile = any(ismember(oldFileNames, tiegcmCurrentFiles(i).name));
    if ~oldsContainThisFile
        newFiles(i) = true;
    end
end

end

function [thisFileDatenums] = giveTiegcmDatenums(tiegcmFile)

modelYear = ncread(tiegcmFile, 'year')'; modelYear = modelYear(1);
modelTimes = double(ncread(tiegcmFile, 'mtime'));
modelDatenums = repmat(datenum(num2str(modelYear),'yyyy'), 1, size(modelTimes,2));
modelDatenums = modelDatenums + modelTimes(1,:) + modelTimes(2,:)/24 + modelTimes(3,:)/1440 - 1;
thisFileDatenums = modelDatenums';

end