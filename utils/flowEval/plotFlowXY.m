function f = plotFlowXY( flow, mag )
%PLOTFLOWXY Plot the x,y components of scene flow separately

f = figure; 
flowx = flow(:, :, 1);
flowy = flow(:, :, 2);
if nargin < 2
    mag = max([max(abs(flowx(:))),max(abs(flowy(:)))]);
    mag = max(mag, eps);
end

subplot(1,2,1); imagesc(flowx,[-mag mag]); colorbar; colormap(hot(256));
pbaspect([size(flowx,2) size(flowy,1) 1]);
title('X Motion');
set(gca,'XTick',[]);
set(gca,'XTickLabel',[]);
set(gca,'YTick',[]);
set(gca,'YTickLabel',[]);

subplot(1,2,2); imagesc(flowy,[-mag mag]); colorbar; colormap(hot(256));
pbaspect([size(flowx,2) size(flowy,1) 1]);
title('Y Motion');
set(gca,'XTick',[]);
set(gca,'XTickLabel',[]);
set(gca,'YTick',[]);
set(gca,'YTickLabel',[]);

set(gcf, 'Position', [100, 100, 1000, 400])
end

