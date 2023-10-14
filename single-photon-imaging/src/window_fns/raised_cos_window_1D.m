function [window] = raised_cos_window_1D(L)
%RAISED_COS_WINDOW Raised cosine window
    window = 0.5 - 0.5 * cos(2 * pi * (double(0:L-1)' + 0.5) / L);
end