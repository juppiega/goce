function tleMap = downloadTLEs(objectIDs, beginDatenums, endDatenums)

load Bfactors.dat

satellites = Bfactors(:,1);
B = Bfactors(:,2);
sig_B = Bfactors(:,3);

if length(beginDatenums) == 1
    beginDatenums = beginDatenums * ones(size(objectIDs));
end

if length(endDatenums) == 1
    endDatenums = endDatenums * ones(size(objectIDs));
end

tleMap = containers.Map('KeyType', 'double', 'ValueType', 'any'); % Sorted by ascending epoch (oldest first).

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
    tleLine1 = fgetl(tleFile);
    if length(tleLine1) ~= 69
        warning(['Something wrong with TLE download or object ',num2str(objectIDs(i)),' does not have any elements between given dates.'])
        fclose(tleFile);
        continue;
    end
    
    tleLine2 = fgetl(tleFile);
    satrec = twoline2rv(72, tleLine1, tleLine2, '0', '0');
    
    [~,output] = system(['wc -l ', outputName]);
    c = strsplit(output);
    numTLEs = str2num(c{1}) / 2;
    
    sgp4SatInfos = repmat(satrec, numTLEs, 1);
    
    for j = 2:numTLEs
        tleLine1 = fgetl(tleFile);
        tleLine2 = fgetl(tleFile);
        sgp4SatInfos(j) = twoline2rv(72, tleLine1, tleLine2, '0', '0');
    end
    
    Bind = satellites == objectIDs(i);
    if sum(Bind) == 0
        error(['Could not find Btrue for object: ', objectIDs(i)]);
    end
    if sum(Bind) > 1
        error(['File contained Btrue twice (or more) for object: ', objectIDs(i), '. A satellite must have only one Btrue value!']);
    end
    Btrue = B(Bind);
    sig_Btrue = sig_B(Bind);
    
    satelliteData = struct('sgp4info', sgp4SatInfos, 'Btrue', Btrue, 'sig_Btrue', sig_Btrue);
    tleMap(objectIDs(i)) = satelliteData;
    
    fclose(tleFile);
    delete(outputName);
end

end