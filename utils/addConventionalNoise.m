function imn = addConventionalNoise(imsum, maxFlux, tau, eta, readNoiseSigma, darkCurrent, fullWell, imgScale, bitDepth, clipNegative)
%ADDCONVENTIONALNOISE Add photon noise and read noise to conventional image
%V2: include dark current noise

if nargin < 10
    clipNegative = true;
end

% Compute the light intensity
imsum = imsum * maxFlux * tau;
if size(imsum, 3) == 1
    eta = mean(eta);
    imsum = imsum * eta;
else
    assert(size(imsum, 3) == 3 && numel(eta) == 3);
    imsum = imsum .* reshape(eta,1,1,3);
end

% Sample
imn = (poissrnd(imsum) + normrnd(0, sqrt(readNoiseSigma^2+darkCurrent*tau), size(imsum)));
if clipNegative
    imn(imn < 0) = 0;
end
imn(imn > fullWell) = fullWell;

if isempty(imgScale)
    imgScale = 1/maxFlux/tau./eta;
end
if size(imn, 3) == 1
    imgScale = mean(imgScale);
    imn = imn * imgScale;
else
    assert(size(imn, 3) == 3 && numel(imgScale) == 3);
    imn = imn .* reshape(imgScale,1,1,3);
end

quanScale = 2 ^ bitDepth;
imn = round(imn * quanScale) / quanScale;
end

