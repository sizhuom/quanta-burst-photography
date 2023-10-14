function [ys] = translate_1D(y, shift, boundary_condition, extrapvals)
%TRANSLATE_1D
    if nargin < 3
        boundary_condition = 'constant';
    end
    
    if nargin < 4
        extrapvals = [0 0];
    end
    
    if isstring(extrapvals) || ischar(extrapvals)
        if strcmpi(extrapvals, 'edge_vals')
            extrapvals = [y(1) y(end)];
        elseif strcmpi(extrapvals, 'edge_mean')
            extrapvals = repmat(0.5 * (y(1) + y(end)), [1 2]);
        elseif strcmpi(extrapvals, 'mean')
            extrapvals = repmat(mean(y, 'all'), [1 2]);
        end
    end
    
    y = squeeze(y);
    L = numel(y);
    padding = ceil(abs(shift)) + 3;
    t = (-padding : L-1+padding)';
    yp = zeros(L + 2*padding, 1);
    yp(padding+1:padding+L) = y;
    if strcmpi(boundary_condition, 'periodic')
        % fill left
        wend = padding;
        while wend >= L
            wstart = wend - L + 1;
            yp(wstart:wend) = y;
            wend = wend - L;            
        end
        yp(1:wend) = y(end-wend+1:end);
        
        % fill right
        wstart = padding + L + 1;
        while wstart <= 2*padding + 1
            wend = wstart + L - 1;
            yp(wstart:wend) = y;
            wstart = wstart + L;
        end
        yp(wstart:end) = y(1:(1+L+2*padding-wstart));
    elseif strcmpi(boundary_condition, 'constant')
        yp(1:padding) = extrapvals(1);
        yp(padding+L+1:end) = extrapvals(2);
    elseif strcmpi(boundary_condition, 'reflect')
        if padding <= L
            yp(padding:-1:1) = y(1:padding);
            yp(end:-1:padding+L+1) = y(end-padding+1:end);
        else
            yp(padding:-1:padding-L+1) = y;
            yp(padding+2*L:-1:padding+L+1) = y;
            yp(1:padding-L) = extrapvals(1);
            yp(padding+2*L+1:end) = extrapvals(2);
        end
    end
    
    ts = t - shift;
    ys = interp1(t, yp, ts, 'spline');
    ys = ys(padding+1:padding+L);
end