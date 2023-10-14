function [shift] = optflow_2D_b(I1, I2, mIx, mIy)
%OPTFLOW_2D_B
    if nargin < 3
        [mIx, sIx] = mag_grad_2D_b(I1, 2);
    end
    
    if nargin < 4
        [mIy, sIy] = mag_grad_2D_b(I1, 1);
    end
    
    mIx(sIx) = -mIx(sIx);
    mIy(sIy) = -mIy(sIy);
    It = xor(I1, I2);
    sIt = I1;
    It(sIt) = -It(sIt);
    shift = -[reshape(mIx, [], 1) reshape(mIy, [], 1)] \ reshape(2 * It, [], 1);
end

