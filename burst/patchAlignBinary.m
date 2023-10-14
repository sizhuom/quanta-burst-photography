function [flows, flowrs] = patchAlignBinary(imbs, param)
%PATCHALIGNBINARY Align binary sequence using patchAlign
% V1: based on patchAlignAggrePyramidV4
%Input:
%  imbs: 1D cell array of binary frames
%  param: struct that contains following fields:
%    alignTWSize: window size for temporal reconstruction
%    alignTWNum: number of temporal windows (determines total number of
%                frames being used), > 2
%    numLevels: number of pyramid levels
%    patchSizes: array that contains patch sizes for each level.
%                Assumption: (1) size of each level is divisible by the
%                patch size, (2) patchSizes(1) must be an even number
%    upsampleRatios: array that contains the upsample ratios for each
%                    pyramid level (upsampled from previous level)
%                    (first entry always 1)
%    searchRadii: array that contains the search radius at each lavel
%    numLKIters: number of Lucas-Kanade iterations for subpixel refinement
%    refFrame: reference frame #
%    imgScale: linear scaling factor for intensity image
%    doRefine: do flow refinement
%    resultDir: directory to save results in
%    debug: whether or not to print debug information
%Output:
%  flows: 1D cell array of computed flows
%  flowrs: 1D cell array of refined flows
ts = tic;
%% Parameters
N = numel(imbs);
[H, W, C] = size(imbs{1});
alignTWSize = param.alignTWSize;
alignTWNum = param.alignTWNum;
if alignTWSize * alignTWNum > N
    error('alignTWSize * alignTNNum must be no greater than N!');
end

refFrame = param.refFrame;
if alignTWSize * alignTWNum < refFrame
    error('alignTWSize * alignTNNum must be no smaller than refFrame!');
end
refBlock = floor((refFrame - 1) / alignTWSize) + 1;
param.refImage = refBlock;

% Get the frame number for block i, frame j (i,j starting from 1)
    function idx = frameIdx(i, j)
        idx = (i - 1) * alignTWSize + j;
    end

imgScale = param.imgScale;

resultDir = param.resultDir;

%% Temporal summing
h = H; w = W;
blockAggres = cell(1, alignTWNum);
for i = 1:alignTWNum
    S = zeros(h, w, param.dataType);
    for j = 1:alignTWSize
        if C == 1
            S = S + imbs{frameIdx(i,j)};
        else
            S = S + mean(imbs{frameIdx(i,j)}, 3);
        end
    end
    blockAggres{i} = S / alignTWSize;
end

if param.debug
    blockRecons = cell(1, alignTWNum);
    for i = 1:alignTWNum
        blockRecons{i} = mleImage(blockAggres{i} * alignTWSize, alignTWSize, imgScale, true);
        imwrite(blockRecons{i}, fullfile(resultDir, sprintf('blockRecons%d-l%d.png', i, 1)));
    end
else
    blockRecons = {};
end

%% Align
flows = patchAlign(blockAggres, param, blockRecons);
if param.doRefine
    flowrs = patchAlignRefine(blockAggres, flows, param, blockRecons);
else
    flowrs = {};
end

fprintf('Alignment done. ');
toc(ts);
end

