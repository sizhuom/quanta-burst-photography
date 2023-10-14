function [B] = anscombeBinomial(A, T)
%ANSCOMBEBINOMIAL Anscombe binomial transform

B = sqrt(T+1/2) * asin(sqrt((A+3/8)/(T+3/4)));

end

