function refFrame = calcRefFrame(twSize, twNum)
%CALCREFFRAME Helper function for computing the center reference frame id 
cenFrame = round((1+twSize)/2);
frameOffset = (1:twSize)-cenFrame;
cenBlock = round((1+twNum)/2);
blockOffset = (1:twNum)-cenBlock;
refFrame = 1-(blockOffset(1)*twSize+frameOffset(1));

end

