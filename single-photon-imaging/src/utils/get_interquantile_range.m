function [iqr] = get_interquantile_range(obs, q1, q2, dim)
%GET_INTERQUANTILE_RANGE
    if nargin < 4
        dim = ndims(obs);
    end
    up_quant = quantile(obs, q2, dim);
    down_quant = quantile(obs, q1, dim);
    iqr = up_quant - down_quant;
end

