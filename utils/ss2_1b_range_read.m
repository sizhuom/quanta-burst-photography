function imbs = ss2_1b_range_read(imgDir, startFrame, endFrame)
%SS2_1B_RANGE_READ Read a continuous subsequence from captured ss2_1b data
%INPUT:
% imgDir: "output_images" folder that contains info.json and part_*.mat
% startFrame, endFrame: 1-based frame idx, inclusive

seqInfo = loadjson(fullfile(imgDir, 'info.json'));
startPart = ceil(startFrame / seqInfo.no_frames);
startOffset = startFrame - (startPart-1) * seqInfo.no_frames;
endPart = ceil(endFrame / seqInfo.no_frames);
endOffset = endFrame - (endPart-1) * seqInfo.no_frames;

imbs = cell(1,(endFrame-startFrame+1));

if startPart == endPart % Only 1 part
    load(fullfile(imgDir, sprintf('part_%d.mat', startPart)), 'OUTPUT');
    for j = startOffset:endOffset
        imbs{j-startOffset+1} = OUTPUT(:,:,j);
    end
else % Multiple parts
    % First part
    load(fullfile(imgDir, sprintf('part_%d.mat', startPart)), 'OUTPUT');
    for j = startOffset:seqInfo.no_frames
        imbs{j-startOffset+1} = OUTPUT(:,:,j);
    end
    
    % Middle parts
    for i = startPart+1:endPart-1
        load(fullfile(imgDir, sprintf('part_%d.mat', i)), 'OUTPUT');
        for j = 1:seqInfo.no_frames
            imbs{(i-startPart)*seqInfo.no_frames-startOffset+1+j} = OUTPUT(:,:,j);
        end
    end
    
    % Final part
    load(fullfile(imgDir, sprintf('part_%d.mat', endPart)), 'OUTPUT');
    for j = 1:endOffset
        imbs{(endPart-startPart)*seqInfo.no_frames-startOffset+1+j} = OUTPUT(:,:,j);
    end

end

