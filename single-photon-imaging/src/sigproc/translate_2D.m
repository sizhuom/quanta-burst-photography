function [Ys] = translate_2D(Y, shift)
%TRANSLATE_2D
    Ys = imtranslate(Y, shift, 'cubic');
end

