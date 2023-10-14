function imd = conventionalBm3d(im, param)
%CONVENTIONALBM3D Apply bm3d to conventional reconstruction result

if param.bm3dSigma == 0
    imd = im;
else
    % Assuming the number of photons can be recovered from parameters
    imgScale = param.imgScale;
    if isempty(imgScale)
        imgScale = 1/param.maxFlux/param.twSize/param.twNum/param.tau./param.eta;
    end
    imgScale = reshape(imgScale, 1, 1, []);
    imNumPhotons = im ./ imgScale;
    % Apply Anscombe binomial transform first
    imt = anscombePoisson(imNumPhotons);
    maxScale = max(imt(:));
    imt = imt / maxScale;
    % Apply bm3d
    if size(imt, 3) > 1
        [~, imtd] = CBM3D(1, imt, param.bm3dSigma*sqrt(255/maxScale), 'np', 0);
    else
        [~, imtd] = BM3D(1, imt, param.bm3dSigma*sqrt(255/maxScale), 'np', 0);
    end
    imd = anscombePoissonInv(imtd*maxScale) .* imgScale;
    imd(imd<0) = 0;
end
end

