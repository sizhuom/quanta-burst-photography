function [imcs, imcSingle] = simConventionalFromIntensity(ims, param, cfa)
%SIMCONVENTIONALFROMINTENSITY Simulate conventional camera images
% from a fast sequence of simulated (noiseless) intensity images

%% Parameters
N = numel(ims);
[H, W, C] = size(ims{1});
if nargin < 3
    cfa = [];
else
    assert(mod(H,size(cfa,1)) == 0 && mod(W,size(cfa,2)) == 0);
    assert(C == 3 && size(cfa,3) == 3);
    cfaFull = repmat(cfa, H/size(cfa,1), W/size(cfa,2), 1);
end
twSize = param.twSize;
twNum = param.twNum;

if twSize * twNum > N
    error('twSize * twNum must be no greater than N!');
end

refFrame = param.refFrame;
if twSize * twNum < refFrame
    error('twSize * twNum must be no smaller than refFrame!');
end

% Get the frame number for block i, frame j (i,j starting from 1)
    function idx = frameIdx(i, j)
        idx = (i - 1) * twSize + j;
    end

maxFlux = param.maxFlux;
tau = param.tau;
eta = param.eta;

%% Simulate captured images
ts = tic;
imcs = cell(1, twNum);
imcSingle = zeros(H, W, C);
for i = 1:twNum
    imtemp = zeros(H, W, C);
    for j = 1:twSize
        imtemp = imtemp + ims{frameIdx(i,j)};
    end
    imtemp = imtemp / twSize;
    imcSingle = imcSingle + imtemp;
    imtemp = addConventionalNoise(imtemp, maxFlux, twSize*tau, eta,...
        param.readNoiseSigma, param.darkCurrent, param.fullWell, param.imgScale, param.bitDepth, param.clipNegative);
    if ~isempty(cfa)
        imtemp = sum(imtemp .* cfaFull, 3);
    end
    imcs{i} = imtemp;
end
imcSingle = imcSingle / twNum;
imcSingle = addConventionalNoise(imcSingle, maxFlux, twSize*twNum*tau, eta,...
    param.readNoiseSigma, param.darkCurrent, param.fullWell, param.imgScale, param.bitDepth, param.clipNegative);
if ~isempty(cfa)
    imcSingle = sum(imcSingle .* cfaFull, 3);
end
toc(ts);

end

