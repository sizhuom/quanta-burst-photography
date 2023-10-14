function [S] = genSumImg(I, T, maxFlux, tau, eta, dcr)
%GENSUMIMG Generate a sum image from an intensity image and imaging
%parameters, which is the sum of T binary images

if nargin < 4
    tau = 1;
end
if nargin < 5
    eta = 1;
end
if nargin < 6
    dcr = 0;
end

P = 1 - exp(-I*maxFlux*tau.*reshape(eta,1,1,[]) - dcr*tau);
S = binornd(T, P);

end

