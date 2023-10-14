function [flows, flowrs] = conventionalPatchAlign(imcs, param)
%CONVENTIONALPATCHALIGN Align conventional images
%Input:
%  imcs: 1D cell array of images, (float, normalized to [0, 1])
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
%    doRefine: do flow refinement
%    resultDir: directory to save results in
%    debug: whether or not to print debug information
%Output:
%  flows: 1D cell array of computed flows
%  imds: 1D cell array of demosaickd images
%  flowrs: 1D cell array of refined flows
ts = tic;

%% Align
flows = patchAlign(imcs, param);
if param.doRefine
    flowrs = patchAlignRefine(imcs, flows, param);
else
    flowrs = {};
end

fprintf('Alignment done. ');
toc(ts);
end