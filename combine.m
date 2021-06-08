% Specify the folder where the files live.
myFolder = 'C:\Users\kolleggerm1\Desktop\HMSPdryZ';
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.dat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
 
  data = load(fullFileName);
  C{k} = data; % store the data
end

A = cat(3,C{:}); % make the 3d matrix


















% % [filename directory_name] = uigetfile('*.dat', 'MultiSelect','on');
% 
% fullname = fullfile('C:\Users\kolleggerm1\Desktop\HMSPdryZ', '*.dat');
% N = 471 ;%numel(S);
% C = cell(1,N);
% for k = 1:N
%      data = load(fullname);
%      C{k} = data; % store the data
% end
% % A = cat(3,C{:}); % concatenate into 3D array