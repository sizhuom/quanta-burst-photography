function [w] = raised_cos_window_2D(M, N)
%RAISED_COS_WINDOW_2D
    wv = raised_cos_window_1D(M);
    wh = raised_cos_window_1D(N);
    w = wv * wh';
end

