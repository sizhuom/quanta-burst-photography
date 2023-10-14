function [C] = genSingleBitJot(I, maxFlux, tau, eta, dcr, readNoiseSigma, thresh)
%GENSINGLEBITJOT Generate single bit jot counts

if nargin < 3
    tau = 1;
end
if nargin < 4
    eta = 1;
end
if nargin < 5
    dcr = 0;
end
if nargin < 6
    readNoiseSigma = 0;
end
if nargin < 7
    thresh = 1;
end

% P = poissrnd(I*maxFlux*tau.*reshape(eta,1,1,[]) + dcr*tau) + normrnd(zeros(size(I)), readNoiseSigma);
P = poissrnd(I*maxFlux*tau.*reshape(eta,1,1,[]) + dcr*tau) - dcr*tau + normrnd(zeros(size(I)), readNoiseSigma);
C = P > thresh-0.5;

end

