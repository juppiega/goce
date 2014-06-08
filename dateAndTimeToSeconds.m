
function [timestampInSeconds] = dateAndTimeToSeconds(dates, hms, dateFormat, timeSeparator)
% dateAndTimeToSeconds(dates, hms, dateFormat, timeSeparator)
% returns dates and HoursMinutesSeconds into one single timestamp

standardDates = datenum(dates, dateFormat);
HMSinSeconds = hms2sec(hms, timeSeparator);

secondsInDay = 60 * 60 * 24;

timestampInSeconds = standardDates * secondsInDay + HMSinSeconds;
timestampInSeconds = timestampInSeconds - timestampInSeconds(1);
end


function [HMSinSeconds] = hms2sec(hms, timeSeparator)
% hms2sec(hms, timeSeparator)

hmsSize = size(hms);
hmsLength = hmsSize(1);
HMSinSeconds = zeros(hmsLength, 1);

parfor i = 1:hmsLength
    
    hmsSplitted = strsplit(hms{i}, timeSeparator);
    
    hours = str2double(hmsSplitted(1));
    minutes = str2double(hmsSplitted(2));
    seconds = str2double(hmsSplitted(3));
    
    HMSinSeconds(i) = 3600 * hours + 60 * minutes + seconds; 
end

end
