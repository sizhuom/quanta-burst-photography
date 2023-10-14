function [Lambda, sigma2] = mleIntensity(S, T, tau, eta, fixInf)
%MLEINTENSITY MLE of intensity

if nargin < 3
    tau = 1;
end
if nargin < 4
    eta = 1;
end
if nargin < 5
    fixInf = true;
end

if fixInf
    if numel(T) > 1
        S(S == T) = T(S == T) - 1;
    else
        S(S == T) = T - 1;
    end
end
Lambda = -log(1 - S ./ T) / tau / eta;

if nargout > 1
    LambdaFixed = Lambda;
    LambdaFixed(LambdaFixed == 0) = -log(1 - 1 ./ T) / tau / eta;
    sigma2 = (exp(LambdaFixed*tau*eta) - 1) / (T*(tau*eta)^2);
end

end

