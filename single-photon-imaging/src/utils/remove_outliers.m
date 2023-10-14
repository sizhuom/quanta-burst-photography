function [Dg, inlier_idxs, Dgstack] = remove_outliers(D, sigma)
%REMOVE_OUTLIERS
    assert(size(D, 1) == size(D, 2));
    have_outliers = true;
    inlier_idxs = 1:size(D, 1);
    Dgstack = {};
    Dg = D;
    Dgstack = [Dgstack Dg];
    while have_outliers
        oidxs = outlier_idxs(Dg, 2, sigma);
        inlier_idxs(oidxs) = [];
        if sum(oidxs) == 0
            have_outliers = false;
        else
            Dg = Dg(~oidxs,:);
            Dg = Dg(:, ~oidxs);
            Dgstack = [Dgstack Dg];
        end
    end
end
