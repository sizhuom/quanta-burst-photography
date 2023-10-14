function [flows] = patchAlign(ims, param, imv)
%PATCHALIGN Align intensity images using patch matching
% V2: only do L1 LK flow when doSR
%Input:
%  ims: 1D cell array of images, (float, normalized to [0, 1])
%  param: struct that contains following fields:
%    refImage: index of the reference image in ims
%    numLevels: number of pyramid levels
%    patchSizes: array that contains patch sizes for each level.
%                Assumption: (1) size of each level is divisible by the
%                patch size, (2) patchSizes(1) must be an even number
%    upsampleRatios: array that contains the upsample ratios for each
%                    pyramid level (upsampled from previous level)
%                    (first entry always 1)
%    searchRadii: array that contains the search radius at each lavel
%    numLKIters: number of Lucas-Kanade iterations for subpixel refinement
%    resultDir: directory to save results in
%    debug: whether or not to print debug information
%  imv: images used for debugging (warped imv will be saved)
%Output:
%  flows: 1D cell array of computed flows
%TODO: adapt to RGB images (lkAlign)

%% Parameters
[H, W, C] = size(ims{1});
N = numel(ims);
refImage = param.refImage;

if nargin < 3 || isempty(imv)
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

resultDir = param.resultDir;

%% Align blocks
numLevels = param.numLevels;
patchSizes = param.patchSizes;
patchStride = patchSizes(1) / 2; % force this to be half patch size
upsampleRatios = param.upsampleRatios;
searchRadii = param.searchRadii;
numStrides = patchSizes(1) / patchStride;
assert(~mod(numStrides, 1));
hs = (H-patchSizes(1))/patchStride+1;
ws = (W-patchSizes(1))/patchStride+1;
assert(~mod(hs, 1) && ~mod(ws, 1));
flows = cell(1, N);

% Pyramid for center block
P0 = buildAggrePyramid(img{refImage}, upsampleRatios);

% Loop over all blocks
for i = 1:N
    if i == refImage
        flows{i} = zeros([hs ws 2]);
        continue
    end
    if param.debug
        fprintf('Block %d: ', i);
    end
    timeBlockStart = tic;
    
    % Build pyramid
    P1 = buildAggrePyramid(img{i}, upsampleRatios);
    
    % Coarse-to-fine matching
    hl = size(P1{numLevels},1) / patchSizes(numLevels);
    wl = size(P1{numLevels},2) / patchSizes(numLevels);
    assert(~mod(hl, 1) && ~mod(wl, 1));
    initMatch = zeros([hl wl 1 2], param.dataType);
    for l = numLevels:-1:2
        if param.debug
            fprintf('L%d:', l);
        end
        bestMatch = zeros([hl wl 2], param.dataType);
        % Find the best match at current level
        for j = 1:hl
            for k = 1:wl
                ylb = 1+(j-1)*patchSizes(l); % y lower bound
                xlb = 1+(k-1)*patchSizes(l); % x lower bound
                bestScore = Inf;
                for m = 1:size(initMatch, 3)
                    [currMatch, currScore] = blockMatch2d(P0{l}, P1{l}, [ylb xlb], patchSizes(l), searchRadii(l), initMatch(j,k,m,:));
                    if currScore < bestScore
                        bestMatch(j,k,:) = currMatch;
                        bestScore = currScore;
                    end
                end
                % Use LK to refine 
                tempf = lkAlign(P0{l}(ylb:ylb+patchSizes(l)-1,xlb:xlb+patchSizes(l)-1), P1{l},...
                    param.numLKIters, bestMatch(j,k,:)+cat(3,xlb-1,ylb-1));
                bestMatch(j,k,:) = tempf - [xlb-1 ylb-1];
            end
        end
        
        % Upsample the match results
        if l > 2
            hl = size(P1{l-1},1) / patchSizes(l-1);
            wl = size(P1{l-1},2) / patchSizes(l-1);
            assert(~mod(hl, 1) && ~mod(wl, 1));
            initMatch = zeros([hl wl 3 2], param.dataType);
            bestMatch = round(bestMatch * upsampleRatios(l));
%             yv = 0:hl-1;
%             xv = 0:wl-1;
%             yr = floor((yv)/upsampleRatios(l)) + 1;
%             xr = floor((xv)/upsampleRatios(l)) + 1;
            yv = (0:hl-1)*patchSizes(l-1)/upsampleRatios(l)/patchSizes(l);
            xv = (0:wl-1)*patchSizes(l-1)/upsampleRatios(l)/patchSizes(l);
            yr = floor(yv) + 1;
            xr = floor(xv) + 1;
            initMatch(:,:,1,1) = bestMatch(yr,xr,1);
            initMatch(:,:,1,2) = bestMatch(yr,xr,2);
            
            yrn = yr - 1;
            masky = mod(yv,1) >= 1/2 & yr < size(bestMatch,1) | yr == 1;
            yrn(masky) = yr(masky) + 1;
            initMatch(:,:,2,1) = bestMatch(yrn,xr,1);
            initMatch(:,:,2,2) = bestMatch(yrn,xr,2);
            
            xrn = xr - 1;
            maskx = mod(xv,1) >= 1/2 & xr < size(bestMatch,2) | xr == 1;
            xrn(maskx) = xr(maskx) + 1;
            initMatch(:,:,3,1) = bestMatch(yr,xrn,1);
            initMatch(:,:,3,2) = bestMatch(yr,xrn,2);
        end
        
        % Debug information
        if param.debug
            flowwarp = repelem(bestMatch,H/size(bestMatch,1),W/size(bestMatch,2),1)/size(P1{l-1},1)*H;
            flowhsv = drawFlowHSV(flowwarp);
            imwrite(flowhsv, fullfile(resultDir, sprintf('flow%d-l%d.png', i, l)));
            fprintf('%g ', toc(timeBlockStart));
        end
    end
    
    % Finest level
    % Upsample
    if param.debug
        fprintf('L1:');
    end
    initMatch = zeros([hs ws 3 2], param.dataType);
    if numLevels > 1
        bestMatch = round(bestMatch * upsampleRatios(2));
        yv = (0:patchStride:H-patchSizes(1))/upsampleRatios(2)/patchSizes(2);
        xv = (0:patchStride:W-patchSizes(1))/upsampleRatios(2)/patchSizes(2);
        yr = floor(yv) + 1;
        xr = floor(xv) + 1;
        initMatch(:,:,1,1) = bestMatch(yr,xr,1);
        initMatch(:,:,1,2) = bestMatch(yr,xr,2);
        
        yrn = yr - 1;
        masky = mod(yv,1) >= 1/2 & yr < size(bestMatch,1) | yr == 1;
        yrn(masky) = yr(masky) + 1;
        initMatch(:,:,2,1) = bestMatch(yrn,xr,1);
        initMatch(:,:,2,2) = bestMatch(yrn,xr,2);
        
        xrn = xr - 1;
        maskx = mod(xv,1) >= 1/2 & xr < size(bestMatch,2) | xr == 1;
        xrn(maskx) = xr(maskx) + 1;
        initMatch(:,:,3,1) = bestMatch(yr,xrn,1);
        initMatch(:,:,3,2) = bestMatch(yr,xrn,2);
    end
    
    % Match
    bestMatch = zeros([hs ws 2], param.dataType);
    % Find the best match at current level
    for j = 1:hs
        for k = 1:ws
            ylb = 1+(j-1)*patchStride; % y lower bound
            xlb = 1+(k-1)*patchStride; % x lower bound
            bestScore = Inf;
            for m = 1:size(initMatch, 3)
                [currMatch, currScore] = blockMatch2d(P0{1}, P1{1}, [ylb xlb], patchSizes(1), searchRadii(1), initMatch(j,k,m,:));
                if currScore < bestScore
                    bestMatch(j,k,:) = currMatch;
                    bestScore = currScore;
                end
            end
            % Use LK to refine
            if ~(param.fastMode && ~param.doSR)
                tempf = lkAlign(P0{1}(ylb:ylb+patchSizes(1)-1,xlb:xlb+patchSizes(1)-1), P1{1},...
                    param.numLKIters, bestMatch(j,k,:)+cat(3,xlb-1,ylb-1));
                bestMatch(j,k,:) = tempf - [xlb-1 ylb-1];
            end
        end
    end
    flows{i} = bestMatch;
    
    % Debug information
    if param.debug
        [X, Y] = meshgrid(1:W, 1:H);
        flowwarp = repelem(bestMatch(1:numStrides:end,1:numStrides:end,:),patchSizes(1),patchSizes(1),1);
        flowhsv = drawFlowHSV(flowwarp);
        imwrite(flowhsv, fullfile(resultDir, sprintf('flow%d-l%d.png', i, 1)));
        imvWarped = zeros(H, W, size(imv{i},3), param.dataType);
        for c = 1:size(imv{i},3)
            imvWarped(:,:,c) = interp2(imv{i}(:,:,c), X+flowwarp(:,:,1), Y+flowwarp(:,:,2), 'cubic');
        end
        imwrite(imvWarped, fullfile(resultDir, sprintf('imWarped%d-l%d.png', i, 1)));
        fprintf('%g\n', toc(timeBlockStart));
    end
end

end
