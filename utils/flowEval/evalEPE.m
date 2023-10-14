function [ epe ] = evalEPE( flow, gtFlow, borderSize )
%EVALEPE Evaluate endpoint error

if nargin < 7
    borderSize = 0;
end
se = (flow(:,:,1)-gtFlow(:,:,1)).^2 + (flow(:,:,2)-gtFlow(:,:,2)).^2;
se = se(1+borderSize:end-borderSize,1+borderSize:end-borderSize);
epe = mean(sqrt(se(:)));

end

