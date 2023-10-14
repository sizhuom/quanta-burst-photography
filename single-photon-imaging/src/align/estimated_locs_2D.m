function [Dhat] = estimated_locs_2D(D)
%ESTIMATED_LOCS_2D
    assert(size(D, 1) == size(D, 2) && size(D, 3) == 2);
    T = size(D, 1);
    Dhat = zeros(T, 2);
    for t = 2:T
        Dhat(t,:) = nan;
        for tp = t-1:-1:1
            if ~any(isnan(D(tp,t,:)))
                if any(isnan(Dhat(tp,:)))
                    Dhat(tp,:) = 0;
                end
                Dhat(t,:) = Dhat(tp,:) + reshape(D(tp, t, :), [1 2]);
            end
        end        
    end
end

