function [B] = genBinImgCFA(I, maxFlux, tau, eta, dcr, cfa)
%GENBINIMGCFA Generate a binary image from an intensity image and imaging
%parameters with color filter array
% cfa: HxWx3 array, HxW is the size of cfa atom (usually H=W), see the
%      Bayer examples (e.g. 'rggb') below

if ischar(cfa)
    switch cfa
        case 'rggb'
            cfa = cat(3, [1 0; 0 0], [0 1; 1 0], [0 0; 0 1]);
        case 'bggr'
            cfa = cat(3, [0 0; 0 1], [0 1; 1 0], [1 0; 0 0]);
        case 'grbg'
            cfa = cat(3, [0 1; 0 0], [1 0; 0 1], [0 0; 1 0]);
        case 'gbrg'
            cfa = cat(3, [0 0; 1 0], [1 0; 0 1], [0 1; 0 0]);
    end
end
assert(mod(size(I,1),size(cfa,1)) == 0 && mod(size(I,2),size(cfa,2)) == 0);
assert(size(I,3) == 3 && size(cfa,3) == 3);
cfaFull = repmat(cfa, size(I,1)/size(cfa,1), size(I,2)/size(cfa,2), 1);
Icfa = sum(I .* reshape(eta,1,1,[]) .* cfaFull, 3);

P = exp(-Icfa*maxFlux*tau - dcr*tau);
B = rand(size(P)) > P;


    
end

