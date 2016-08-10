function S = computeDst(S)

dstFile = fopen('DST.txt');
if dstFile == -1
    error('dstFile open unsuccesful')
end

dstData = textscan(dstFile, '%s %s %f %f', 'MultipleDelimsAsOne',1);

try
    timestampsDstDatenum = datenum(strcat(dstData{1}, dstData{2}), 'yyyy-mm-ddHH:MM:SS.FFF');
catch
    fclose(dstFile);
    system('sed -i.bak ''/|/d'' ./DST.txt');
    dstFile = fopen('DST.txt');
    dstData = textscan(dstFile, '%s %s %f %f', 'MultipleDelimsAsOne',1);
    timestampsDstDatenum = datenum(strcat(dstData{1}, dstData{2}), 'yyyy-mm-ddHH:MM:SS.FFF');
end
dst = dstData{4};

S.dst = interp1(timestampsDstDatenum, dst, S.timestamps);

end