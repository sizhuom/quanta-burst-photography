function [imr] = postMerge(S, param, isSR, dcr)
%POSTMERGE Post-merge processing
%V3: 
% 201230: fix sigma

if nargin < 3
    isSR = false;
end
if isSR
    T = param.srTWNum * param.srTWSize;
else
    T = param.mergeTWNum * param.mergeTWSize;
end

% Remove dark counts
if nargin > 3 && ~isempty(dcr)
    if param.correctDCR
        S = S - dcr*param.tau*T;
    end
end
S(S<0) = 0;

% Denoise
if param.bm3dSigma > 0
% Apply Anscombe binomial transform first
    imt = anscombeBinomial(S, T);
    maxScale = max(imt(:));
    imt = imt / maxScale;
    if size(S, 3) > 1
        [~, imtd] = CBM3D(1, imt, param.bm3dSigma/maxScale, 'np', 0);
    else
        [~, imtd] = BM3D(1, imt, param.bm3dSigma/maxScale, 'np', 0);
    end
    imd = anscombeBinomialInv(imtd*maxScale, T);
else
    imd = S;
end

% Invert the response curve
imr = mleImage(imd, T, param.imgScale, true); % fix inf here
imr = real(imr);
imr(imr < 0) = 0;
