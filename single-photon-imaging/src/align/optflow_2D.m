function [shift] = optflow_2D(I1, I2, sigma, grad_threshold)
%OPTFLOW_2D
    if nargin < 3
        sigma = 0;
    end
    
    if nargin < 4
        grad_threshold = 1e-4;
    end
    
    if ~isfloat(I1)
        I1 = double(I1);
    end
    
    if ~isfloat(I2)
        I2 = double(I2);
    end
    
    if sigma ~= 0
        I1 = imgaussfilt(I1, sigma);
        I2 = imgaussfilt(I2, sigma);
    end
    
    I1_nextx = circshift(I1, -1, 2);
%     I1_prevx = circshift(I1, 1, 2);
    I1_nexty = circshift(I1, -1, 1);
%     I1_prevy = circshift(I1, 1, 1);
%     Ix = (I1_nextx - I1_prevx) / 2;
%     Iy = (I1_nexty - I1_prevy) / 2;
    Ix = I1_nextx - I1;
    Iy = I1_nexty - I1;
    It = I2 - I1;
    
    Ix = reshape(Ix(1:end-1, 1:end-1), [], 1);
    Iy = reshape(Iy(1:end-1, 1:end-1), [], 1);
    It = reshape(It(1:end-1, 1:end-1), [], 1);
%     Ix = reshape(diff(I1(1:end-1,:), 1, 2), [], 1);
%     Iy = reshape(diff(I1(:,1:end-1), 1, 1), [], 1);
%     It = reshape(I2(1:end-1, 1:end-1) - I1(1:end-1, 1:end-1), [], 1);
%     valid_idxs = Ix .^ 2 + Iy .^2 >= grad_threshold^2;
%     if sum(valid_idxs) > 0
%         Ix = Ix(valid_idxs);
%         Iy = Iy(valid_idxs);
%         It = It(valid_idxs);
    shift = -[Ix, Iy] \ It;
%     else
%         shift = [0 0];
%     end
end

