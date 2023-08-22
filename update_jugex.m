updateInstalledVersion();


% Update the installed version of export_fig from the latest version online
function updateInstalledVersion()
    % Download the latest version of export_fig into the export_fig folder
    zipFileName = 'https://github.com/INM-1-FZJ-Cytoarchitecture/jugex/archive/refs/heads/main.zip';
    fprintf('Downloading latest version of %s from %s...\n', mfilename, zipFileName);
    folderName = fileparts(which(mfilename('fullpath')));
    targetFileName = fullfile(folderName, datestr(now,'yyyy-mm-dd.zip'));
    try
        folder = hyperlink(['matlab:winopen(''' folderName ''')'], folderName);
    catch  % hyperlink.m is not properly installed
        folder = folderName;
    end
    try
        urlwrite(zipFileName,targetFileName); %#ok<URLWR>
    catch err
        error('export_fig:update:download','Error downloading %s into %s: %s\n',zipFileName,targetFileName,err.message);
    end

    % Unzip the downloaded zip file in the export_fig folder
    fprintf('Extracting %s...\n', targetFileName);
    try
        unzip(targetFileName,folderName);
        % Fix issue #302 - zip file uses an internal folder export_fig-master
        subFolder = fullfile(folderName,'export_fig-master');
        try movefile(fullfile(subFolder,'*.*'),folderName, 'f'); catch, end %All OSes
        try movefile(fullfile(subFolder,'*'),  folderName, 'f'); catch, end %MacOS/Unix
        try movefile(fullfile(subFolder,'.*'), folderName, 'f'); catch, end %MacOS/Unix
        try rmdir(subFolder); catch, end
    catch err
        error('export_fig:update:unzip','Error unzipping %s: %s\n',targetFileName,err.message);
    end

    % Notify the user and rehash
    fprintf('Successfully installed the latest %s version in %s\n', mfilename, folder);
    clear functions %#ok<CLFUNC>
    rehash
end