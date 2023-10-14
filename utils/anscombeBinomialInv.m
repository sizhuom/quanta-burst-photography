function [A] = anscombeBinomialInv(B, T)
%ANSCOMBEBINOMIALINV Inverse Anscombe binomial transform

A = (sin(B/sqrt(T+1/2))).^2 * (T+3/4) - 3/8;

end

