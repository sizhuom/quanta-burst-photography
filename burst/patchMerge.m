function [imo] = patchMerge(patches, param)
%PATCHMERGE Merge pre-aligned intensity patches via patchwise Wiener filtering
%V1: based on patchWienerMerge V5
%Input:
%  patches: 4D array of pre-aligned blocks of size
%       (hs*patchSize, ws*patchSize, C, N)
%  param: struct that contains following fields:
%    refImage: index of the reference image in ims
%    patchSizes: array that contains patch sizes for each level.
%    wienerC: tuning parameter C for wiener filtering
%    debug: whether or not to print debug information
%Output:
%  S: Sum image

%% Parameters
[Hp, Wp, C, N] = size(patches);
H = param.H;
W = param.W;
patchSize = param.patchSizes(1);
patchStride = patchSize / 2; % force this to be half patch size

refImage = param.refImage;

winWeights = raised_cos_window_2D(patchSize, patchSize);

%% Merge color channels sequentially
imo = zeros(H, W, C, param.dataType);
hs = Hp/patchSize;
ws = Wp/patchSize;

for c = 1:C
    % Swap the central block and the first block
    tmp = patches(:,:,c,1);
    patches(:,:,c,1) = patches(:,:,c,refImage);
    patches(:,:,c,refImage) = tmp;
    
    % Then merge all blocks
    Sa = zeros(H, W, param.dataType);
    accWeights = zeros(H, W, param.dataType);
    fprintf('Block merging...\n');
    for i = 1:hs
        if param.debug
            fprintf('%d', i);
        end
        for j = 1:ws
            ylb = 1+(i-1)*patchSize; % y lower bound
            xlb = 1+(j-1)*patchSize; % x lower bound
            Sdenoised = wiener_denoise_t(squeeze(patches(ylb:ylb+patchSize-1,xlb:xlb+patchSize-1,c,:)), param.wienerC);
            ylb = 1+(i-1)*patchStride; % y lower bound
            xlb = 1+(j-1)*patchStride; % x lower bound
            Sa(ylb:ylb+patchSize-1,xlb:xlb+patchSize-1) = Sa(ylb:ylb+patchSize-1,xlb:xlb+patchSize-1)+Sdenoised.*winWeights;
            accWeights(ylb:ylb+patchSize-1,xlb:xlb+patchSize-1) = accWeights(ylb:ylb+patchSize-1,xlb:xlb+patchSize-1) + winWeights;
            if param.debug
                fprintf('.');
            end
        end
        if param.debug
            fprintf('\n');
        end
    end
    Sa = Sa ./ accWeights;
    imo(:,:,c) = Sa;
end
imo(imo<0) = 0;

end

