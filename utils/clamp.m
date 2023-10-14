function [A] = clamp(A, minVal, maxVal)
%CLAMP Clamp values in an array to a given range

if nargin > 1 && ~isempty(minVal)
    A(A<minVal) = minVal;
end
if nargin > 2 && ~isempty(maxVal)
    A(A>maxVal) = maxVal;
end
end

