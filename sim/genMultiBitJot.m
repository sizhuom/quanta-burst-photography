function [C] = genMultiBitJot(I, numBits, maxFlux, tau, eta, dcr, readNoiseSigma)
%GENMULTIBITJOT Generate multi-bit jot counts

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

% P = poissrnd(I*maxFlux*tau.*reshape(eta,1,1,[]) + dcr*tau) + normrnd(zeros(size(I)), readNoiseSigma);
P = poissrnd(I*maxFlux*tau.*reshape(eta,1,1,[]) + dcr*tau) - dcr*tau + normrnd(zeros(size(I)), readNoiseSigma);
C = round(P);
L = 2^numBits - 1;
C = C / L;
C(C<0) = 0;
C(C>1) = 1;

end

