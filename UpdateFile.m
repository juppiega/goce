function UpdateFile(filename)
%UpdateFile Checks for the named file and updates it if necessary.
%   This function is designed to download the text files necessary to run
%   JB2008. It can only run with MATLAB 2014b or later.
%
% UpdateFile('filename.txt') saves the filename in the current directorty
% if the file is missing or out of date. There are no output arguments.



if exist(filename,'file') == 2
    fileInfo = dir(filename);
    if floor(fileInfo.datenum) > datenum(date)
        websave(filename,strcat('http://sol.spacenvironment.net/jb2008/indices/',filename));
        fprintf('Downloading the latest version of %s\n',filename)
    end
elseif exist(filename,'file') ~= 2
    websave(filename,strcat('http://sol.spacenvironment.net/jb2008/indices/',filename));
    fprintf('Downloading the latest version of %s\n',filename)
end








end

