function [mgI, sgI] = mag_grad_2D_b(I, dim, diff_type)
%MAG_GRAD_2D_b
% Sign is given by I itself. 1 => negative gradient, 0 => positive gradient
    if nargin < 3
        diff_type = 'central';
    end
    
    I_next = circshift(I, -1, dim);
    
    if strcmpi(diff_type, 'central')
        I_prev = circshift(I, 1, dim);
        mgI = xor(I_prev, I_next); % actually half of this value
        if nargout == 2
            sgI = I_prev;
        end
    elseif strcmpi(diff_type, 'forward')
        mgI = xor(I, I_next);
        if nargout == 2
            sgI = I;
        end
    end
end