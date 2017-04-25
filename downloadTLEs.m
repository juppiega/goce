function tleMap = downloadTLEs(objectIDs, beginDatenums, endDatenums, ignoreBtrue)
pause on

downloadFreq = 19; % per minute
dt = 60 / downloadFreq / 86400;
persistent previousTime
if isempty(previousTime)
    previousTime = 0;
end

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
    t = now;
    
    beginString = datestr(beginDatenums(i), 'yyyy-mm-dd-HH:MM:SS');
    endString = datestr(endDatenums(i), 'yyyy-mm-dd-HH:MM:SS');
    
    if (t - previousTime) < dt
        pause((previousTime + dt - t)*86400);
    end
    command = ['python tleDownloader.py ',num2str(objectIDs(i)),' ',...
                beginString,' ',endString, ' tleAccess.txt'];
    disp(objectIDs(i))
    previousTime = now;
      
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
    if (nargin <= 3 || ~ignoreBtrue) && sum(Bind) == 0
        error(['Could not find Btrue for object: ', nu2str(objectIDs(i))]);
    end
    if sum(Bind) > 1
        error(['File contained Btrue twice (or more) for object: ', objectIDs(i), '. A satellite must have only one Btrue value!']);
    end
    if nargin <= 3 || ~ignoreBtrue
        Btrue = B(Bind);
        sig_Btrue = sig_B(Bind);
    else
        Btrue = 0.1;
        sig_Btrue = 0.01;
    end
    
    satelliteData = struct('sgp4info', sgp4SatInfos, 'Btrue', Btrue, 'sig_Btrue', sig_Btrue);
    tleMap(objectIDs(i)) = satelliteData;
    
    fclose(tleFile);
    delete(outputName);
end

end