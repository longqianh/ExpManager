% rename 'tif' files to 'mat'

% Set the directory path where your files are located
directory = 'E:\Longqian\Experiments\20230517\out_17500_48x4\scan_out';

% Create a file pattern to match the files you want to rename
filePattern = '*.tif'; % Modify this pattern to match your file extension or name

% Get a list of all files that match the file pattern in the directory
fileList = dir(fullfile(directory, filePattern));

% Iterate through each file and modify the file name
for i = 1:numel(fileList)
    % Get the current file name
    currentName = fileList(i).name;
    tmp=split(currentName,'.');
    % Modify the file name as desired
    newName = [tmp{1},'.mat']; % Modify this line to change the new file name format
    
    % Create the full paths for the old and new file names
    oldFullPath = fullfile(directory, currentName);
    newFullPath = fullfile(directory, newName);
    
    % Rename the file
    movefile(oldFullPath, newFullPath);
end
