function [ima, S] = naiveRecons(imbs, param)
%NAIVERECONS Naive reconstruction by simple summing the binary images and
%compute MLE
%V2: use mleImage instead of mleIntensity, ignores absolute intensity
%Input:
%  imbs: 1D cell array of binary frames
%  param: struct that contains following fields:
%    mergeTWSize: window size for temporal reconstruction
%    mergeTWNum: number of temporal windows (determines total number of
%                frames being used), > 2
%    refFrame: reference frame #
%    imgScale: linear scaling factor for intensity image
%    debug: whether or not to print debug information
%Output:
%  ima: reconstructed image
%  S: sum image

%% Parameters
N = numel(imbs);
[H, W, C] = size(imbs{1});
twSize = param.mergeTWSize;
twNum = param.mergeTWNum;

if twSize * twNum > N
    error('twSize * twNum must be no greater than N!');
end

% Get the frame number for block i, frame j (i,j starting from 1)
    function idx = frameIdx(i, j)
        idx = (i - 1) * twSize + j;
    end

imgScale = param.imgScale;

%% Sum and compute
S = zeros(H, W, C, param.dataType);
for i = 1:twNum
    for j = 1:twSize
        S = S + imbs{frameIdx(i,j)};
    end
end
ima = mleImage(S, twNum*twSize, imgScale);

end
