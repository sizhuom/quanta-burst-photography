function elem = randelem(arr, sz)
%RANDELEM Randomly pick an element from an array

if nargin < 2
    sz = 1;
end
elem = arr(randi(numel(arr), sz));

end

