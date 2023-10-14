function [psnr] = evalPSNR(imest,imgt)
%EVALPSNR Evaluate PSNR
% 201217: fix a bug: use mean(_, 'all') instead of mean2()

assert(isfloat(imest) && isfloat(imgt));
imgt(imgt<0) = 0; imgt(imgt>1) = 1;
imest(imest<0) = 0; imest(imest>1) = 1;
mse = mean((imgt - imest).^2, 'all');
psnr = 10*log10(1/mse);

end

