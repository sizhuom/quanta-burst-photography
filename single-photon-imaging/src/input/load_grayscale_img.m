function [img] = load_grayscale_img(img_path)
    img = imread(img_path);
    if size(img, 3) > 1
        img = rgb2gray(img);
    end
    img = im2double(img);
end

