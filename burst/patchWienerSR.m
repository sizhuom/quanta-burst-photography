function S = patchWienerSR(imbs, flows, param, mergeSumImage)
%PATCHWIENERSR Super-resolution reconstruction from binary images
% v6: minor fixes, clean up the code a little bit
% Change parameter to mergeSumImage.
%Input:
%  imbs: 1D cell array of binary frames
%  flows: patchwise flow fields for keyframes
%  param: struct that contains following fields:
%    patchSize: spatial patch size for wiener filtering (divisible)
%    alignTWSize: temporal window size for alignment
%    alignTWNum: number of temporal windows for alignment
%    mergeTWSize: temporal window size for merging
%    mergeTWNum: number of temporal windows for merging
%    srTWSize: temporal window size for SR
%    srTWNum: number of temporal windows for SR, total number of
%                frames used for SR <= for alignment/merging
%    refFrame: reference frame #
%    srScale: output scale, outputSize = ceil(originalSize*srScale)
%    combineRadius: radius of the neighborhood to look up for merging
%    imgScale: linear scaling factor for the recovered intensity image
%    k_detail, k_denoise, k_stretch, k_shrink, D_th, D_tr: kernel
%                                                          parameters
%    wienerSRC: tuning parameter C for wiener filtering
%    debug: whether or not to print debug information
%  mergeSumImage: sum image generated by patchWienerMerge 
%Output:
%  S: sum image

ts0 = tic;
%% Parameters
[H, W, C] = size(imbs{1});
srScale = param.srScale;
Ho = ceil(H*srScale);
Wo = ceil(W*srScale);
patchSize = param.patchSizes(1);
patchStride = patchSize / 2; % force this to be half patch size
alignTWSize = param.alignTWSize;
alignTWNum = param.alignTWNum;
srTWSize = param.srTWSize;
srTWNum = param.srTWNum;

if srTWSize * srTWNum > alignTWSize * alignTWNum
    error('srTWSize * srTWNum must be no greater than alignTWSize * alignTNNum!');
end

refFrame = param.refFrame;
if srTWSize * srTWNum < refFrame
    error('srTWSize * srTWNum must be no smaller than refFrame!');
end

% "reference frame" in each sr block
cenFrame = mod(refFrame - 1, srTWSize) + 1; 
refBlock = floor((refFrame - 1) / srTWSize) + 1;

% Get the frame number for block i, frame j (i,j starting from 1)
    function idx = frameIdx(i, j)
        idx = (i - 1) * srTWSize + j;
    end

% Get the align block and frame subscripts for a give nframe no.
    function [i, j] = alignSub(idx)
        i = floor((idx-1)/alignTWSize)+1;
        j = mod(idx-1, alignTWSize) + 1;
    end

imgScale = param.imgScale;

if nargin < 4
    mergeSumImage = [];
end


%% Merge
% "reference frame" in each align block
alignCenFrame = mod(refFrame - 1, alignTWSize) + 1; 

% Preprocess the flows for interpolation
flowsr = cell(1, alignTWNum+2);
for i = 1:alignTWNum
    flowsr{i+1} = flows{i};
end
flowsr{1} = 2*flowsr{2} - flowsr{3};
flowsr{alignTWNum+2} = 2*flowsr{alignTWNum+1} - flowsr{alignTWNum};

% Function that returns the interpolated flow for a frame
    function iflow = interpFlow(i, j)
        idx = frameIdx(i, j);
        [ai, aj] = alignSub(idx);
        assert(ai > 0 && ai <= alignTWNum && aj > 0 && aj <= alignTWSize);
        if aj < alignCenFrame
            iflow = (alignCenFrame-aj)/alignTWSize*flowsr{ai} + (aj+alignTWSize-alignCenFrame)/alignTWSize*flowsr{ai+1};
        else
            iflow = (alignCenFrame+alignTWSize-aj)/alignTWSize*flowsr{ai+1} + (aj-alignCenFrame)/alignTWSize*flowsr{ai+2};
        end
    end

% Phase 1: similar to patchWienerMerge
fprintf('Phase 1...\n');
hs = (H-patchSize)/patchStride+1;
ws = (W-patchSize)/patchStride+1;
xv = repelem((0:ws-1)*patchStride,1,patchSize)+repmat(1:patchSize,1,ws);
yv = repelem((0:hs-1)*patchStride,1,patchSize)+repmat(1:patchSize,1,hs);
[X, Y] = meshgrid(xv, yv);
Sr = zeros(H, W, C, param.dataType);
blockMerge = cell(1, C);
warpX = zeros([hs*patchSize, ws*patchSize, srTWNum], param.dataType);
warpY = zeros([hs*patchSize, ws*patchSize, srTWNum], param.dataType);
winWeights = raised_cos_window_2D(patchSize, patchSize);
weightMap = repmat(winWeights, hs, ws);
if param.debug
    ybv = repelem((0:H/patchSize-1)*2*patchSize,1,patchSize)+repmat(1:patchSize,1,H/patchSize);
    xbv = repelem((0:W/patchSize-1)*2*patchSize,1,patchSize)+repmat(1:patchSize,1,W/patchSize);
end
for c = 1:C
    ts = tic;
    % First merge each block into a single image
    % (This is an over-complete representation)
    blockMerge{c} = zeros(hs*patchSize, ws*patchSize, srTWNum, param.dataType);
    fprintf('Averaging within each block (Channel %d)...\n', c);
    for i = 1:srTWNum
        Sb = zeros(hs*patchSize, ws*patchSize, param.dataType);
        countMap = zeros(size(Sb), param.dataType);
        if param.debug
            fprintf('%d', i);
        end
        for j = 1:srTWSize
            curFlow = interpFlow(i, j);
            flowwarp = repelem(curFlow, patchSize, patchSize, 1);
            destX = X + round(flowwarp(:,:,1));
            destY = Y + round(flowwarp(:,:,2));
            imbwarped = interp2(double(imbs{frameIdx(i,j)}(:,:,c)), destX, destY, 'nearest');
            countMap = countMap + isfinite(imbwarped);
            imbwarped(~isfinite(imbwarped)) = 0;
            Sb = Sb + imbwarped;
            
            if c == 1 && j == cenFrame
                warpX(:,:,i) = destX - flowwarp(:,:,1);
                warpY(:,:,i) = destY - flowwarp(:,:,2);
            end
            if param.debug
                fprintf('.');
            end
        end
        if param.debug
            fprintf('\n');
        end
        countMap(countMap == 0) = 1;
        Sb = Sb ./ countMap;
        blockMerge{c}(:,:,i) = Sb;
    end
    % Swap the central block and the first block
    % Use the merge image instead of the first block
    tmp = blockMerge{c}(:,:,1);
    if ~isempty(mergeSumImage)
        blockMerge{c}(:,:,1) = interp2(mergeSumImage(:,:,c)/param.mergeTWSize/param.mergeTWNum, X, Y, 'nearest');
    else
        blockMerge{c}(:,:,1) = blockMerge{c}(:,:,refBlock);
    end
    blockMerge{c}(:,:,refBlock) = tmp;
    if c == 1
        tmp = warpX(:,:,1);
        warpX(:,:,1) = warpX(:,:,refBlock);
        warpX(:,:,refBlock) = tmp;
        tmp = warpY(:,:,1);
        warpY(:,:,1) = warpY(:,:,refBlock);
        warpY(:,:,refBlock) = tmp;
    end
    toc(ts);
    
    % Apply Wiener filter to each block
    ts = tic;
    Sa = zeros(H, W, param.dataType);
    accWeights = zeros(H, W, param.dataType);
    fprintf('Wiener filtering (Channel %d)...\n', c);
    for i = 1:hs
        if param.debug
            fprintf('%d', i);
        end
        for j = 1:ws
            ylb = 1+(i-1)*patchSize; % y lower bound
            xlb = 1+(j-1)*patchSize; % x lower bound
            blockMerge{c}(ylb:ylb+patchSize-1,xlb:xlb+patchSize-1,:) = wiener_denoise_tf(blockMerge{c}(ylb:ylb+patchSize-1,xlb:xlb+patchSize-1,:), param.wienerSRC);
            Sdenoised = mean(blockMerge{c}(ylb:ylb+patchSize-1,xlb:xlb+patchSize-1,:), 3);
            ylb = 1+(i-1)*patchStride; % y lower bound
            xlb = 1+(j-1)*patchStride; % x lower bound
            Sa(ylb:ylb+patchSize-1,xlb:xlb+patchSize-1) = Sa(ylb:ylb+patchSize-1,xlb:xlb+patchSize-1)+Sdenoised .* winWeights;
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
    Sr(:,:,c) = Sa;
    toc(ts);
end
Sr(Sr < 0) = 0;
% if param.debug
%     for c = 1:C
%         for i = 1:srTWNum
%             bm = blockMerge{c}(ybv,xbv,c);
%             bmr = mleImage(bm*srTWSize, srTWSize, imgScale);
%             imwrite(bmr, fullfile(param.resultDir, sprintf('srBlockMergeRecons%d-c%d.png', i, c)));
%         end
%     end
% end

if ~isempty(mergeSumImage)
    imguide = mean(mergeSumImage, 3);
else
    imguide = mean(Sr, 3);
end

% Then merging using anisotropic Gaussian kernel
S = zeros(Ho, Wo, C, param.dataType);
for c = 1:C    
    ts = tic;
    % Compute the structure tensors
    [Ix, Iy] = imgradientxy(imguide, 'sobel');
    Ix = Ix/8; Iy = Iy/8; % normalize
    Ixx = Ix .^ 2; Ixy = Ix .* Iy; Iyy = Iy .^ 2;
    STa = imfilter(Ixx, ones(3));
    STb = imfilter(Ixy, ones(3));
    STd = imfilter(Iyy, ones(3));
    STT = STa + STd; STD = STa .* STd - STb .^ 2;
    lambda1 = STT/2 + sqrt(STT.^2/4-STD);
    lambda2 = STT - lambda1;
    A = min(sqrt(lambda1./lambda2), 5);
    A(isnan(A)) = 1;
    D = 1 - sqrt(lambda1)/param.D_tr + param.D_th;
    D(D<0) = 0; D(D>1) = 1;
    k1hat = param.k_detail * param.k_stretch * A;
    k2hat = param.k_detail ./ (param.k_shrink * A);
    k1 = ((1-D).*k1hat + D*param.k_detail*param.k_denoise).^2;
    k2 = ((1-D).*k2hat + D*param.k_detail*param.k_denoise).^2;
    
    fprintf('Merging (Channel %d)...\n', c);
    % Mapping assumption: 
    % (1) same image size: (0.5, H+0.5) for both images
    % (2) same image center: (1+H)/2 maps to (1+Ho)/2
    iOffset = (1+H)/2 - (1+Ho)/2/srScale;
    jOffset = (1+W)/2 - (1+Wo)/2/srScale;
    for i = 1:Ho
        if param.debug
            fprintf('%d', i);
        end
        for j = 1:Wo
            accC = 0;
            accW = 0;
            iOrig = iOffset + i/srScale; % pixel index in original grid
            jOrig = jOffset + j/srScale;
            iOrigR = round(iOrig);
            jOrigR = round(jOrig);
            % Compute the covariance matrix Omega
            if STb(iOrigR,jOrigR) == 0
                if STa(iOrigR,jOrigR) > STd(iOrigR,jOrigR)
                    e1 = [1; 0];
                    e2 = [0; 1];
                else
                    e1 = [0; 1];
                    e2 = [1; 0];
                end
            else
                e1 = [STb(iOrigR,jOrigR); lambda1(iOrigR,jOrigR)-STa(iOrigR,jOrigR)];
                e1 = e1 / norm(e1);
                e2 = [STb(iOrigR,jOrigR); lambda2(iOrigR,jOrigR)-STa(iOrigR,jOrigR)];
                e2 = e2 / norm(e2);
            end
%             Omega = [e1 e2]*[k1 0;0 k2]*[e1 e2]';
            OmegaInv = [e2 e1]*[1/k1(iOrigR,jOrigR) 0; 0 1/k2(iOrigR,jOrigR)]*[e2 e1]';
            for k = max(1,iOrigR-param.combineRadius):min(H,iOrigR+param.combineRadius)
                for l = max(1,jOrigR-param.combineRadius):min(W,jOrigR+param.combineRadius)
                    % Find the positions of (k,l) in the over-complete
                    % tiles representation
                    % TODO: FOR NOW ASSUME patchSize = 2 * patchStride
                    ky = [floor((k-1)/patchStride)*patchSize+mod(k-1,patchStride)+1,...
                        floor((k-1-patchStride)/patchStride)*patchSize+mod(k-1,patchStride)+1+patchStride]; 
                    lx = [floor((l-1)/patchStride)*patchSize+mod(l-1,patchStride)+1,...
                        floor((l-1-patchStride)/patchStride)*patchSize+mod(l-1,patchStride)+1+patchStride]; 
                    if k <= patchStride
                        ky = ky(1);
                    elseif k > H - patchStride
                        ky = ky(2);
                    end
                    if l <= patchStride
                        lx = lx(1);
                    elseif l > W - patchStride
                        lx = lx(2);
                    end
                    for m = ky
                        for n = lx
                            d = [reshape(warpX(m,n,:)-jOrig, 1, []); reshape(warpY(m,n,:)-iOrig, 1, [])];
                            w0 = exp(-sum(OmegaInv*d.*d,1)/2);
                            w = w0 * weightMap(m,n);
                            w(isnan(w)) = 0;
                            accC = accC + sum(w.*reshape(blockMerge{c}(m,n,:),1,[]));
                            accW = accW + sum(w);
                        end
                    end
                end
            end
            S(i,j,c) = max(accC / accW, 0);
            if param.debug
                fprintf('.');
            end
        end
        if param.debug
            fprintf('\n');
        end
    end
    toc(ts);
end
S = S*srTWNum*srTWSize;
fprintf('Super-resolution done. ');
toc(ts0);

end

