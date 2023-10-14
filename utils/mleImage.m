function [Lambda, sigma2] = mleImage(S, T, scale, fixInf)
%MLEIMAGE Convert a sum image into an linear intensity image

if nargin < 3
    scale = 1;
end
if nargin < 4
    fixInf = false;
end

if fixInf
    if numel(T) > 1
        S(S > T - 1) = T(S > T - 1) - 1;
    else
        S(S > T - 1) = T - 1;
    end
end

Lambda = -log(1 - S ./ T);
if nargin > 2
    if size(Lambda, 3) == 1
        scale = mean(scale);
        Lambda = Lambda * scale;
    else
        assert(size(Lambda, 3) == 3);
        if numel(scale) == 3
            Lambda = Lambda .* reshape(scale,1,1,3);
        else
            Lambda = Lambda * scale;
        end
    end
else
    Lambda = Lambda / max(Lambda(:));
end

if nargout > 1
    LambdaFixed = Lambda;
    LambdaFixed(LambdaFixed == 0) = -log(1 - 1 ./ T) / tau / eta;
    sigma2 = (exp(LambdaFixed*tau*eta) - 1) / (T*(tau*eta)^2);
end

end

