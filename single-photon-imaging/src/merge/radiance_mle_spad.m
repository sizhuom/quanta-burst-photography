function [mle] = radiance_mle_spad(mean_obs, num_obs, total_exposure_time, saturation_value)
%RADIANCE_MLE
    if nargin < 2
        num_obs = 1;
    end
    
    if nargin < 3
        total_exposure_time = num_obs;
    end
    
    if nargin < 4
        saturation_value = log(num_obs) + log(2);
    end
    
    mle = min(-log(1 - mean_obs), saturation_value) / (total_exposure_time / num_obs);    
end

