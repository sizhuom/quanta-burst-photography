function [B] = genBinImg(I, maxFlux, tau, eta, dcr)
%GENBINIMG Generate a binary image from an intensity image and imaging
%parameters

if nargin < 3
    tau = 1;
end
if nargin < 4
    eta = 1;
end
if nargin < 5
    dcr = 0;
end

P = exp(-I*maxFlux*tau.*reshape(eta,1,1,[]) - dcr*tau);
B = rand(size(P)) > P;
    
end

