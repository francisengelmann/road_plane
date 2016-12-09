function [boxPix, boxIndex] = getImBox(im, x, y, boxWidth, boxHeight, onlyInd)

if nargin < 5
    boxHeight = boxWidth;
end

if nargin < 6
    onlyInd = false;
end

imheight = size(im, 1);
imwidth = size(im, 2);

startx = max(floor(x - boxWidth/2), 1);
starty = max(floor(y - boxHeight/2), 1);
endx = min(startx+boxWidth-1, imwidth);
endy = min(starty+boxHeight-1, imheight);
if endx == imwidth || endy == imheight
    %boxPix = imresize(boxPix, [boxSize boxSize]);
    startx = endx - boxWidth + 1;
    starty = endy - boxHeight + 1;
end

boxIndex = {starty:endy, startx:endx};
if ~onlyInd
    % normalize to 0-1
    boxPix = double(im(boxIndex{:}, :)) ./ 255.0;
else
    boxPix = [];
end

end
