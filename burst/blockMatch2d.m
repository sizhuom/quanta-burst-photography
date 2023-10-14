function [bestMatch, bestScore] = blockMatch2d(im0, im1, blockUL, blockSize, searchRadius, initMatch)
%BLOCKMATCH2D Brute-force block match in 2D

if nargin < 6
    initMatch = zeros(1, 2);
end
[h, w] = size(im0);
ylb = blockUL(1);
xlb = blockUL(2);
refPatch = im0(ylb:ylb+blockSize-1,xlb:xlb+blockSize-1,:);
bestScore = inf;
bestMatch = [nan nan];
for v = max(1-ylb,initMatch(2)-searchRadius):min(h-blockSize+1-ylb,initMatch(2)+searchRadius)
    for u = max(1-xlb,initMatch(1)-searchRadius):min(w-blockSize+1-xlb,initMatch(1)+searchRadius)
        curPatch = im1(ylb+v:ylb+v+blockSize-1,xlb+u:xlb+u+blockSize-1,:);
        curScore = sum(abs(curPatch-refPatch), 'all');
        if curScore < bestScore
            bestScore = curScore;
            bestMatch = [u v];
        end
    end
end
end

