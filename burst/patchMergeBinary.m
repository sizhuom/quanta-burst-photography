function [S] = patchMergeBinary(imbs, flows, param)
%PATCHMERGEBINARY Merging binary sequence via patchMerge
%V2: 
%  210818: Implement warpTWSize: do not warp every binary frame
%Input:
%  imbs: 1D cell array of binary frames
%  flows: coarse flows from alignment algorithm
%  param: struct that contains following fields:
%    patchSizes: array that contains patch sizes for each level.
%                Assumption: (1) size of each level is divisible by the
%                patch size, (2) patchSizes(1) must be an even number
%                Note: only patchSizes(1) is used in merge
%    alignTWSize: temporal window size for alignment
%    alignTWNum: number of temporal windows for alignment
%    mergeTWSize: temporal window size for merging
%    mergeTWNum: number of temporal windows for merging, total number of
%                frames used for merging <= for alignment
%    refFrame: reference frame #
%    imgScale: linear scaling factor for the recovered intensity image
%    wienerC: tuning parameter C for wiener filtering
%    resultDir: directory to save results in
%    debug: whether or not to print debug information
%Output:
%  S: merged sum image
ts0 = tic;
%% Parameters
[H, W, C] = size(imbs{1});
patchSize = param.patchSizes(1);
patchStride = patchSize / 2; % force this to be half patch size
alignTWSize = param.alignTWSize;
alignTWNum = param.alignTWNum;
mergeTWSize = param.mergeTWSize;
mergeTWNum = param.mergeTWNum;

if mergeTWSize * mergeTWNum > alignTWSize * alignTWNum
    error('mergeTWSize * mergeTWNum must be no greater than alignTWSize * alignTNNum!');
end

refFrame = param.refFrame;
if mergeTWSize * mergeTWNum < refFrame
    error('mergeTWSize * mergeTWNum must be no smaller than refFrame!');
end

refBlock = floor((refFrame - 1) / mergeTWSize) + 1;
param.refImage = refBlock;

% Get the frame number for block i, frame j (i,j starting from 1)
    function idx = frameIdx(i, j)
        idx = (i - 1) * mergeTWSize + j;
    end

% Get the align block and frame subscripts for a give nframe no.
    function [i, j] = alignSub(idx)
        i = floor((idx-1)/alignTWSize)+1;
        j = mod(idx-1, alignTWSize) + 1;
    end

imgScale = param.imgScale;

%% Merge
if alignTWNum == 1
    assert(mergeTWNum == 1 && alignTWSize == mergeTWSize);
    S = zeros(size(imbs{1}), param.dataType);
    for i = 1:alignTWSize
        S = S + imbs{frameIdx(1,i)};
    end
    fprintf('Merging done. ');
    toc(ts0);
    return
end

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

% First merge each block into a single image
ts = tic;
hs = (H-patchSize)/patchStride+1;
ws = (W-patchSize)/patchStride+1;
xv = repelem((0:ws-1)*patchStride,1,patchSize)+repmat(1:patchSize,1,ws);
yv = repelem((0:hs-1)*patchStride,1,patchSize)+repmat(1:patchSize,1,hs);
[X, Y] = meshgrid(xv, yv);
if param.debug
    ybv = repelem((0:H/patchSize-1)*2*patchSize,1,patchSize)+repmat(1:patchSize,1,H/patchSize);
    xbv = repelem((0:W/patchSize-1)*2*patchSize,1,patchSize)+repmat(1:patchSize,1,W/patchSize);
end
blockPatches = zeros(hs*patchSize, ws*patchSize, C, mergeTWNum, param.dataType);
for c = 1:C
    fprintf('Pre-warping Channel %d...\n', c);
    for i = 1:mergeTWNum
        Sb = zeros(hs*patchSize, ws*patchSize, param.dataType);
        countMap = zeros(size(Sb), param.dataType);
        if param.debug
            fprintf('.');
        end
        if isfield(param, 'warpTWSize') && param.warpTWSize > 1
            warpTWSize = param.warpTWSize;
            assert(~mod(mergeTWSize, warpTWSize));
            warpTWNum = mergeTWSize / warpTWSize;
            for j = 1:warpTWNum
                curFlow = interpFlow(i, (j-1)*warpTWSize+round((warpTWSize+1)/2));
                flowwarp = repelem(curFlow, patchSize, patchSize, 1);
                curFrame = zeros(H, W, param.dataType);
                for k = 1:warpTWSize
                    curFrame = curFrame + imbs{frameIdx(i, (j-1)*warpTWSize+k)}(:,:,c);
                end
                imbwarped = interp2(curFrame, X+flowwarp(:,:,1), Y+flowwarp(:,:,2), 'linear');
                countMap = countMap + isfinite(imbwarped) * warpTWSize;
                imbwarped(~isfinite(imbwarped)) = 0;
                Sb = Sb + imbwarped;
            end
        else
            for j = 1:mergeTWSize
                curFlow = interpFlow(i, j);
                flowwarp = repelem(curFlow, patchSize, patchSize, 1);
                curFrame = double(imbs{frameIdx(i,j)}(:,:,c));
                if param.fastMode
                    imbwarped = interp2(curFrame, X+flowwarp(:,:,1), Y+flowwarp(:,:,2), 'nearest');
                else
                    imbwarped = interp2(curFrame, X+flowwarp(:,:,1), Y+flowwarp(:,:,2), 'linear');
                end
                countMap = countMap + isfinite(imbwarped);
                imbwarped(~isfinite(imbwarped)) = 0;
                Sb = Sb + imbwarped;
            end
        end
        countMap(countMap == 0) = 1;
        Sb = Sb ./ countMap;
        
        blockPatches(:,:,c,i) = Sb;
        if param.debug
%             bm = Sb(ybv,xbv);
%             bmr = mleImage(bm*mergeTWSize, countMap(ybv,xbv), imgScale);
%             imwrite(bmr, fullfile(param.resultDir, sprintf('blockMergeRecons%d-c%d.png', i, c)));
        end
    end
    if param.debug
        fprintf('\n');
    end
end
toc(ts);

%% Then Wiener Merge
param.H = H; param.W = W;
S = patchMerge(blockPatches, param);
S = S * mergeTWNum * mergeTWSize;
fprintf('Merging done. ');
toc(ts0);
end

