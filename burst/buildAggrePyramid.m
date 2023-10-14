function P = buildAggrePyramid(ims, upsampleRatios)
%BUILDAGGREPYRAMID Build aggregate image pyramid. Coarser levels are built
%by aggregating photons from finer levels.
%V2: assume images are already normalized to 0...1

numLevels = numel(upsampleRatios);
P = cell(1, numLevels);
P{1} = ims;
S0 = ims;
aggreRatio = 1;
for i = 2:numLevels
    scaleFilter = ones(upsampleRatios(i));
    scaleCenter = floor((upsampleRatios(i)+1)/2);
    S0 = imfilter(S0, scaleFilter);
    S0 = S0(scaleCenter:upsampleRatios(i):end,scaleCenter:upsampleRatios(i):end,:);
    aggreRatio = aggreRatio * upsampleRatios(i)^2;
    P{i} = S0 / aggreRatio;
end

