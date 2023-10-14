function [C] = genSingleBitJotCFA(I, cfa, maxFlux, tau, eta, dcr, readNoiseSigma, thresh)
%GENSINGLEBITJOTCFA Generate single bit jot counts with CFA

if nargin < 4
    tau = 1;
end
if nargin < 5
    eta = 1;
end
if nargin < 6
    dcr = 0;
end
if nargin < 7
    readNoiseSigma = 0;
end
if nargin < 8
    thresh = 1;
end

Icfa = sum(I .* reshape(eta,1,1,[]) .* cfa, 3);
P = poissrnd(Icfa*maxFlux*tau + dcr*tau) - dcr*tau + normrnd(zeros(size(Icfa)), readNoiseSigma);
C = P > thresh-0.5;

end

