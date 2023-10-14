function [binary_img] = get_spad_shot(influx, exposure_time)
%GET_SPAD_SHOT Generate a binary (SPAD model) image from an incoming light flux.
    if nargin < 2
        exposure_time = 1;
    end
    binary_img = false(size(influx));
    binary_img(rand(size(influx)) > exp(-influx * exposure_time)) = true;
end