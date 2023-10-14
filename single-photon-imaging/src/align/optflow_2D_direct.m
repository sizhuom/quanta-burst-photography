function [shift] = optflow_2D_direct(Ix, Iy, It)
%OPTFLOW_2D_DIRECT
    ix = reshape(Ix, [], 1);
    iy = reshape(Iy, [], 1);
    it = reshape(It, [], 1);
    keep = ix ~= 0 | iy ~= 0;
    if numel(find(keep, 1)) > 0
        shift = -[ix(keep) iy(keep)] \ it(keep);
    else
        shift = [0; 0];
    end
end