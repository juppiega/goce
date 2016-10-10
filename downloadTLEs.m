function tleMap = downloadTLEs(objectIDs, beginDatenums, endDatenums)

if length(beginDatenums) == 1
    beginDatenums = beginDatenums * ones(size(objectIDs));
end

if length(endDatenums) == 1
    endDatenums = endDatenums * ones(size(objectIDs));
end

for i = 1:length(objectIDs)
    beginString = datestr(beginDatenums(i), 'yyyy-mm-dd-HH:MM:SS');
    endString = datestr(endDatenums(i), 'yyyy-mm-dd-HH:MM:SS');
    
    command = ['python3 tleDownloader.py ',num2str(objectIDs(i)),' ',...
                beginString,' ',endString, ' tleAccess.txt'];
            
    [status, output] = system(command);
    if status ~= 0
        error(['Something wrong with python script tleDownloader.py. System output: ', output])
    end
    
    beginString = datestr(beginDatenums(i), 'yyyy-mm-dd-HHMMSS');
    endString = datestr(endDatenums(i), 'yyyy-mm-dd-HHMMSS');
    
    outputName = ['TLE_',sprintf('%05d',objectIDs(i)),'_',beginString,'--',endString,'.txt'];
    tleFile = fopen(outputName);
    tleLine = fgetl(tleFile);
    if length(tleLine) ~= 69
        warning(['Object ',num2str(objectIDs(i)),' does not have any elements between given dates.'])
    end
    
    satrec = twoline2rv( 72, ...
                       longstr1, longstr2, typerun, typeinput);
end

end