function [imcm] = conventionalPatchMerge(imcs, flows, param)
%CONVENTIONALPATCHMERGE Merge conventional camera images using
% This function warp the images and divide them into overlapping patches
% and pass them to patchMerge
%Input:
%  imcs: 1D cell array of images, (float, normalized to [0, 1])
%  param: struct that contains following fields:
%    refImage: index of the reference image in ims
%    patchSizes: array that contains patch sizes for each level.
%    wienerC: tuning parameter C for wiener filtering
%    debug: whether or not to print debug information
%Output:
%  S: Sum image
ts0 = tic;
%% Parameters
[H, W, C] = size(imcs{1});
patchSize = param.patchSizes(1);
patchStride = patchSize / 2; % force this to be half patch size
N = numel(imcs);

%% Merge
ts = tic;
hs = (H-patchSize)/patchStride+1;
ws = (W-patchSize)/patchStride+1;
xv = repelem((0:ws-1)*patchStride,1,patchSize)+repmat(1:patchSize,1,ws);
yv = repelem((0:hs-1)*patchStride,1,patchSize)+repmat(1:patchSize,1,hs);
[X, Y] = meshgrid(xv, yv);

% Warp the demosaicked images
blockPatches = zeros(hs*patchSize, ws*patchSize, C, N);
for c = 1:C
    for i = 1:N
        curFlow = flows{i};
        flowwarp = repelem(curFlow, patchSize, patchSize, 1);
        bsWarped = interp2(imcs{i}(:,:,c), X+flowwarp(:,:,1), Y+flowwarp(:,:,2), 'linear');
        bsWarped(~isfinite(bsWarped)) = 0;
        
        blockPatches(:,:,c,i) = bsWarped;
    end
end
toc(ts);

%% Then Wiener Merge
param.H = H; param.W = W;
imcm = patchMerge(blockPatches, param);
fprintf('Merging done. ');
toc(ts0);
end

