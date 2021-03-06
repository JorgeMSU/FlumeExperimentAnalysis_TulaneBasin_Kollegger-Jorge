% This script is not required to reproduce the results included in Kolleger
% et al. The 3D matrices "Anew" or "LMLP" generated by this script
% are also provided in the Zenodo code repository.   
myFolder = 'C:\topo';
% The "topo" folder should include the files from the SEAD repository (see
% link in acknowledgements or Supplementary Information)associated with:
% hours 50-540 in the case LMLP portion of the experiment
% hours 680-1170 in the case HMSP portion of the experiment
n= 439;% The domain size is nxn in both the LMLP and the HMSP scenarios.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.dat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
m=length(theFiles);
C=zeros(n,n,m);
for k = 1 : m
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
 
  data = load(fullFileName);
  C(:,:,k) = data; % store the data
end
LMLP =C; %Select when appropriate
%Anew =C; %Select when appropriate
 

