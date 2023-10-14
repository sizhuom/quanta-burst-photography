function initQBP
%INITSPS Initialize the path settings
fprintf('Setting up the paths and global variables for Quanta Burst Photography...');
p = mfilename('fullpath');
[root, ~, ~] = fileparts(p);
addpath(root);
subdirs = {'sim', 'utils', 'burst', 'single-photon-imaging'};
for i = 1:length(subdirs)
    addpath(genpath(fullfile(root,subdirs{i})));
end

% JSONLab: https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab-a-toolbox-to-encode-decode-json-files
if isunix
    addpath('~/Documents/MATLAB/fangq-jsonlab-c3eb021');
else
    addpath('D:/MATLAB/fangq-jsonlab-c3eb021');
end

% OpenEXR: https://github.com/skycaptain/openexr-matlab
% https://github.com/edgarv/hdritools
if isunix
    addpath('~/Documents/MATLAB/openexr-matlab');
else
    addpath('D:/lib/HDRITools-0.5.0/matlab');
end

% BM3D V2.0: http://www.cs.tut.fi/~foi/GCF-BM3D/
if isunix
    addpath('~/Documents/MATLAB/BM3D');
else
    addpath('D:/MATLAB/BM3D');
end

fprintf('Done.\n');

end

