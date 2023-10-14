function [Lambda] = mleImageMultiBitJot(S, T, scale)
%MLEIMAGEMULTIBITJOT Convert a sum image into an linear intensity image


Lambda = S / T;
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


end

