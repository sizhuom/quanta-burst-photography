function [shift] = phasecorr_dft_1D(Y1, Y2, method, ignore_higher_freqs)
%PHASECORR_DFT_1D
    if nargin < 3
        method = 'plane_fit';
    end
    
    if nargin < 4
        ignore_higher_freqs = false;
    end
    
    assert(numel(Y1) == numel(Y2));
    L = numel(Y1);
    
    if strcmpi(method, 'plane_fit')
        Y1(1) = 0;
        Y2(1) = 0;
    end
    
    if ignore_higher_freqs
        Y1(3:end-1) = 0;
        Y2(3:end-1) = 0;
    end
    
    CP = Y2 .* conj(Y1);
    Mag_CP = abs(CP);
    t1 = L; t2 = (L*L) * 10; % depends on dynamic range of signal
    valid_freqs = and(Mag_CP >= t1, Mag_CP <= t2);
    Norm_CP = CP ./ abs(CP);
    Norm_CP(~valid_freqs) = 0;
    
    if strcmpi(method, 'plane_fit')
        f = fftshift((0:L-1)' - L/2) * (1 / L) * 2 * pi;
        phase_shifts = angle(Norm_CP(valid_freqs));
        f = f(valid_freqs);
        shift = -(f \ phase_shifts);
    elseif strcmpi(method, 'inv_ft')
        inv_cp = fftshift(idft_1D(Norm_CP, @box_window_1D));
        [m, i_max] = max(inv_cp);
        shift = i_max - 1 - L/2;
        right = 1 + mod(i_max, L);
        left = 1 + mod(i_max - 2, L);
        % Find each shift independently then add it up.
        m_l = inv_cp(left);
        m_r = inv_cp(right);
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
        shift = shift + delta_l + delta_r;
    end
end

