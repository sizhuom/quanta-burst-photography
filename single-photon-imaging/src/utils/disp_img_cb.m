function h = disp_img_cb(Y)
%DISP_IMG_CB
    h0 = figure;
    image(Y, 'CDataMapping', 'scaled');
    colormap gray;
    colorbar;
    if nargout > 0
        h = h0;
    end
end
