
function [TDTdata] = fp2mat(fppath, filename)

% Changes TDT format to *.mat and saves the filename.mat in same folder:
%  INPUTS
%    fppath: folder where TDT data is stored
%    filename: select a name (including .mat) for the file
%  OUTPUT
%    datafile: loads TDTdata to current workspace


TDTdata = TDTbin2mat(fppath);

TDTdata.path = fppath;

% Use the same path to save data.mat:

filename = strcat(fppath,'/',filename);
save(filename,'TDTdata');

% Load TDTdata in workspace:

load(filename, 'TDTdata')

end

