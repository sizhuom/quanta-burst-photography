function [dcr] = genDcrPattern(shape, dcData, scaleMedian)
%GENDCRPATTERN Generate dark count rate fixed pattern
% from measured dark counts
% scaleMedian: scale the overall distribution by scaleMedian /
% median(dcData)

dcData = dcData(:) / (130001 / 97700);
dcData = sort(dcData);
if nargin > 2
    dcData = dcData / median(dcData) * scaleMedian;
end
N = numel(dcData);

ind = randi([1, N], shape);
dcr = dcData(ind);

end

