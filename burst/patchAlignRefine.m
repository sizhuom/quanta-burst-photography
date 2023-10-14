function [flowrs] = patchAlignRefine(ims, flows, param, imv)
%PATCHALIGNREFINE Refine patch alignment
%Based on patchAlignFlowRefineV3
%Input:
%  ims: 1D cell array of frames, (float, normalized to [0, 1])
%  flows: 1D cell array of computed flows
%  param: struct that contains following fields:
%    refImage: index of the reference image in ims
%    patchSizes: array that contains patch sizes for each level.
%                Assumption: (1) size of each level is divisible by the
%                patch size, (2) patchSizes(1) must be an even number
%    flowLambda: weighting factor for the smoothness term, larger lambda =
%                smoother flow field
%    resultDir: directory to save results in
%    debug: whether or not to print debug information
%  imv: images used for debugging (warped imv will be saved)
%Output:
%  flowrs: 1D cell array of refined flows

%% Parameters
[H, W, C] = size(ims{1});
N = numel(ims);
refImage = param.refImage;
resultDir = param.resultDir;

if nargin < 4
    imv = ims;
end
if C == 1
    img = ims;
else
    img = cell(1, N);
    for i = 1:N
        img{i} = rgb2gray(ims{i});
    end
end

%% Align blocks
ts = tic;
patchSizes = param.patchSizes;
patchStride = patchSizes(1) / 2; % force this to be half patch size
numStrides = patchSizes(1) / patchStride;
assert(~mod(numStrides, 1));
hs = (H-patchSizes(1))/patchStride+1;
ws = (W-patchSizes(1))/patchStride+1;
assert(~mod(hs, 1) && ~mod(ws, 1));
flowrs = cell(1, N);

% Loop over all blocks
for i = 1:N
    if i == refImage
        flowrs{i} = zeros([hs ws 2]);
        continue
    end
    
    % Refine the match using patch flow
    pfParam = {'lambda', param.flowLambda, 'patchSize', patchSizes(1), 'patchStride', patchStride, 'display', true};
    temp = estimate_flow_interface(img{refImage}, img{i}, 'classic+nl-patch', pfParam, flows{i});
    flowrs{i} = temp;
    
    % Debug information
    if param.debug
        [X, Y] = meshgrid(1:W, 1:H);
        flowwarp = repelem(temp(1:numStrides:end,1:numStrides:end,:),patchSizes(1),patchSizes(1),1);
        flowhsv = drawFlowHSV(flowwarp);
        imwrite(flowhsv, fullfile(resultDir, sprintf('refined-flow%d-l%d.png', i, 1)));

        imvWarped = zeros(H, W, size(imv{i},3));
        for c = 1:size(imv{i},3)
            imvWarped(:,:,c) = interp2(imv{i}(:,:,c), X+flowwarp(:,:,1), Y+flowwarp(:,:,2), 'cubic');
        end
        imwrite(imvWarped, fullfile(resultDir, sprintf('refined-imWarped%d-l%d.png', i, 1)));
   
        fprintf('.\n');
    end
end

toc(ts);

end
