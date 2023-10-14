function [is_majority] = majority(Y, dim)
%MAJORITY
    assert(islogical(Y));
    is_majority = sum(Y, dim) > (size(Y, dim) / 2);
end