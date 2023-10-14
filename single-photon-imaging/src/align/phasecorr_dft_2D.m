function [shift] = phasecorr_dft_2D(Y1, Y2, method)
%PHASECORR_DFT_2D
    if nargin < 3
        method = 'plane_fit';
    end
    
    [H1, W1] = size(Y1);
    [H2, W2] = size(Y2);
    assert(H1 == H2 && W1 == W2);
    H = H1; W = W1;
    
    if strcmpi(method, 'plane_fit')
        Y1(1, 1) = 0;
        Y2(1, 1) = 0;
    end
    
    CP = Y2 .* conj(Y1);
    Mag_CP = abs(CP);
    t1 = H * W; t2 = (H^2) *(W^2) * 10;
    valid_freqs = and(Mag_CP >= t1, Mag_CP <= t2);
    Norm_CP = CP ./ Mag_CP;
    Norm_CP(~valid_freqs) = 0;
    
    if strcmpi(method, 'plane_fit')
        fx = fftshift((0:W-1)' - W/2) * (1 / W) * 2 * pi;
        fy = fftshift((0:H-1)' - H/2) * (1 / H) * 2 * pi;
        [Fx, Fy] = meshgrid(fx, fy);
        all_shifts = angle(Norm_CP(valid_freqs));
        Fx = Fx(valid_freqs); Fy = Fy(valid_freqs);
        shift = -[reshape(Fx, [], 1), reshape(Fy, [], 1)] \ reshape(all_shifts, [], 1);
    elseif strcmpi(method, 'inv_ft')
        inv_cp = fftshift(idft_2D(Norm_CP, @box_window_2D));
        [m, ind_max] = max(inv_cp(:));
        [i_max, j_max] = ind2sub(size(inv_cp), ind_max);
        shift(1) = j_max - 1 - W/2;
        shift(2) = i_max - 1 - H/2;
        right = 1 + mod(j_max, W);
        left = 1 + mod(j_max - 2, W);
        up = 1 + mod(i_max - 2, H);
        down = 1 + mod(i_max, H);
        % First shift in x
        % Find each shift independently then add it up.
        m_l = inv_cp(i_max, left);
        m_r = inv_cp(i_max, right);
        if m_l < m && m_r < m
            delta_l = -m_l / (m_l + sign(m_l)*m);
            delta_r = m_r / (m_r + sign(m_r)*m);
        elseif m_l == m && m_r == m
            delta_l = 0;
            delta_r = 0;
        elseif m_l < m
            delta_l = 0;
            delta_r = 0.5;
        elseif m_r < m
            delta_l = -0.5;
            delta_r = 0;
        end        
        shift(1) = shift(1) + delta_l + delta_r;
        % Then shift in y
        m_u = inv_cp(up, j_max);
        m_d = inv_cp(down, j_max);
        if m_u < m && m_d < m
            delta_u = -m_u / (m_u + sign(m_u)*m);
            delta_d = m_d / (m_d + sign(m_d)*m);
        elseif m_u == m && m_d == m
            delta_u = 0;
            delta_d = 0;
        elseif m_u < m
            delta_u = 0;
            delta_d = 0.5;
        elseif m_d < m
            delta_u = -0.5;
            delta_d = 0;
        end        
        shift(2) = shift(2) + delta_u + delta_d;
    end
end
