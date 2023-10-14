function [result] = qbpAccumDemoMono(imbs, param, dcr, imgt)
%QBPACCUMDEMOMONO QBP accumulation of frames demo
% based on qbpPipelineMono

result = struct();
resultDir = param.resultDir;

%% Naive reconstruction with simple averaging
ima = naiveRecons(imbs, param);
if param.imgAutoScale
    [ima, param.imgScale] = autoScaleIntensity(ima, 97);
end
% asParam = param;
% asParam.mergeTWNum = 1;
% refBlock = floor((param.refFrame - 1) / param.mergeTWSize);
% imas = naiveRecons(imbs(refBlock*param.mergeTWSize+1:(refBlock+1)*param.mergeTWSize), asParam);

assert(param.removeHP);
if param.removeHP
    tic;
    for i = 1:numel(imbs)
        imbs{i} = removeHotPixels(imbs{i}, dcr, param.hpThresh); 
    end
    toc;
    fprintf('Finished hot pixel correction.\n');
    
    imaf = naiveRecons(imbs, param);
%     imasf = naiveRecons(imbs(refBlock*param.mergeTWSize+1:(refBlock+1)*param.mergeTWSize), asParam);
    
    refBlock = floor((param.refFrame - 1) / param.mergeTWSize);
    assert(refBlock == 0);
    imafList = cell(1, param.mergeTWNum);
    asParam = param;
    for i = 1:param.mergeTWNum
        asParam.mergeTWNum = i;
        imafList{i} = naiveRecons(imbs, asParam);
    end
else
    imaf = [];
%     imasf = [];
end
fprintf('Finished naive reconstruction.\n');

%% Align
[flows, flowrs] = patchAlignBinary(imbs, param);
if param.debug
    save(fullfile(resultDir, 'patchAlign.mat'), 'flows', 'flowrs', 'param');
end
result.flows = flows;
result.flowrs = flowrs;
fprintf('Finished alignment.\n');

%% Merge
imrList = cell(1, param.mergeTWNum);
for i = 1:param.mergeTWNum
    paramNoBM3D = param;
    paramNoBM3D.bm3dSigma = 0;
    paramNoBM3D.alignTWNum = i;
    paramNoBM3D.mergeTWNum = i;
    Sr = patchMergeBinary(imbs, flows, paramNoBM3D);
    imr = postMerge(Sr, paramNoBM3D, false);
    imrList{i} = imr;
end

if param.imgAutoScale
    [imr, imgScaleCorrect] = autoScaleIntensity(imr, 97);
    param.imgScale = param.imgScale * imgScaleCorrect;
    ima = ima * imgScaleCorrect;
    imaf = imaf * imgScaleCorrect;
    for i = 1:param.mergeTWNum
        imafList{i} = imafList{i} * imgScaleCorrect;
        imrList{i} = imrList{i} * imgScaleCorrect;
    end
end

if param.bm3dSigma > 0
    imrbm = postMerge(Sr, param, false);
    if param.saveImages
        imwrite(lin2rgb(imrbm), fullfile(resultDir, 'bm3d_g22.png'));
    end
else
    imrbm = [];
end


if param.saveImages
%     imwrite(lin2rgb(imr), fullfile(resultDir, 'patchMerge_g22.png'));
%     imwrite(lin2rgb(ima), fullfile(resultDir, 'averageRecons_g22.png'));
    if ~isfield(param, 'removeHP') || param.removeHP
%         imwrite(lin2rgb(imaf), fullfile(resultDir, 'averageReconsHPFixed_g22.png'));
    end
    
    for i = 1:param.mergeTWNum
        imwrite(lin2rgb(imafList{i}), fullfile(resultDir, sprintf('naive_%d.png',i)));
        imwrite(lin2rgb(imrList{i}), fullfile(resultDir, sprintf('qbp_%d.png',i)));
    end
end
if param.debug
    save(fullfile(resultDir, 'patchMerge.mat'), 'param', 'Sr', 'imr', 'imrbm', 'imrList');
    save(fullfile(resultDir, 'naiveRecons.mat'), 'ima', 'imaf', 'imafList');
end
result.ima = ima;
% result.imas = imas;
result.imaf = imaf;
% result.imasf = imasf;
result.imafList = imafList;
result.Sr = Sr;
result.imr = imr;
result.imrbm = imrbm;
result.imrList = imrList;
result.param = param;
fprintf('Finished merging.\n');

%% Refine flow and merge
% if param.doRefine
%     Srr = patchMergeBinary(imbs, flowrs, param);
%     imrr = postMerge(Srr, param, false);
%     if param.debug
%         save(fullfile(resultDir, 'patchMerge_refinedFlow.mat'), 'imrr', 'Srr');
%     end
%     if param.saveImages
%         imwrite(lin2rgb(imrr), fullfile(resultDir, 'patchMerge_refinedFlow_g22.png'));
%     end
%     result.Srr = Srr;
%     result.imrr = imrr;
%     fprintf('Finished flow refinement.\n');
% end
% 
% %% Superresolution
% if param.doSR
%     Ssr = patchWienerSR(imbs, flows, param, Sr);
%     imsr = postMerge(Ssr, param, true);
%     if param.debug
%         save(fullfile(resultDir, 'patchWienerSR.mat'), 'imsr', 'Ssr');
%     end
%     if param.saveImages
%         imwrite(lin2rgb(imsr), fullfile(resultDir, 'patchWienerSR_g22.png'));
%     end
%     result.Ssr = Ssr;
%     result.imsr = imsr;
%     fprintf('Finished super-resolution.\n');
% end
% 
% if param.doRefineSR
%     Ssrr = patchWienerSR(imbs, flowrs, param, Sr);
%     imsrr = postMerge(Ssrr, param, true);
%     if param.debug
%         save(fullfile(resultDir, 'patchWienerSR_refinedFlow.mat'), 'imsrr', 'Ssrr');
%     end
%     if param.saveImages
%         imwrite(lin2rgb(imsrr), fullfile(resultDir, 'patchWienerSR_refinedFlow_g22.png'));
%     end
%     result.Ssrr = Ssrr;
%     result.imsrr = imsrr;
%     fprintf('Finished super-resolution with flow refinement.\n');
% end
% 
% %% Compute PSNR
% if param.computePSNR
%     psnr = struct();
%     [H, W, ~] = size(imr);
%     imgtr = imresize(imgt, [H W]);
%     imgtr(imgtr>1) = 1;
%     imgtr(imgtr<0) = 0;
%     psnr.quantaNaivePSNR = evalPSNR(ima, imgtr);
%     psnr.quantaBurstPSNR = evalPSNR(imr, imgtr);
%     
%     if param.doRefine
%         psnr.quantaBurstRefinedPSNR = evalPSNR(imrr, imgtr);
%     end
%     
%     if param.doSR
%         [Hsr, Wsr, ~] = size(imsr);
%         imgtsr = imresize(imgt, [Hsr Wsr]);
%         psnr.quantaBurstSrPSNR = evalPSNR(imsr, imgtsr);
%     end
%     
%     if param.doRefineSR
%         [Hsr, Wsr, ~] = size(imsr);
%         imgtsr = imresize(imgt, [Hsr Wsr]);
%         psnr.quantaBurstSrRefinedPSNR = evalPSNR(imsrr, imgtsr);
%     end
%     
%     if param.debug
%         savejson('', psnr, fullfile(resultDir, 'quanta.json'));
%     end
%     result.psnr = psnr;
% end
end

