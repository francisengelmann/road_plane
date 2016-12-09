function line = ransacLine(X, Y, Z, maxDist, pk, P)

if nargin < 6
    % Desired probability
    P = 0.99999;
end

if nargin < 5
    % Chance to pick both pixels on the line
    pk = 0.2;
end

if nargin < 4
    % Max 10 cm from line
    maxDist = 0.1;
end

numPix = numel(X);
numSamples = round(log(1-P) / log(1-pk));
% 2 points for every sample
POINTS_PER_SAMP = 2;
replacement = numSamples * POINTS_PER_SAMP > numPix;
samples = randsample(numPix, numSamples * POINTS_PER_SAMP, replacement);
maxNumInliers = -1;
maxInliers = [];
bestAvgDist = inf;
% Return a x line by default
point1 = [0 0 0];
point2 = [1 0 0];
v = 0;
% Get the 3d pos of each pixel in a numPixel x 3 array
% where the columns are (X,Y,Z)
pos3d = cat(2, X(:), Y(:), Z(:));
t1 = 0; t2 = 0; t3 = 0;
for i = 1:numSamples
    tic;
    offset = (i-1) * POINTS_PER_SAMP;
    % Get the pixel number of each sample
    s1 = samples(offset+1);
    s2 = samples(offset+2);
    % Get the position of the pixel
    p1 = pos3d(s1, :);
    p2 = pos3d(s2, :);
    % Get the vector from the line formed by the points
    v1 = p2 - p1;
    mag = norm(v1);
    v1 = v1 ./ mag;
    dists = pointLineDist(pos3d, v1, p1);
    t1 = t1 + toc;
    tic;
    inliers = dists <= maxDist;
    inDists = dists(inliers);
    avgDist = mean(inDists);
    %sqDist = sum(inDists .^ 2);
    numInliers = sum(inliers(:));
    t2 = t2 + toc;
    tic;
    if numInliers > maxNumInliers ...
            || (numInliers == maxNumInliers && avgDist < bestAvgDist)
        % If the inliers are the same and the distance is reduced, then
        % it is a better line and use it
        maxNumInliers = numInliers;
        maxInliers = inliers;
        bestAvgDist = avgDist;
        point1 = p1;
        point2 = p2;
        v = v1;
    end
    t3 = t3 + toc;
end
inliers = maxInliers;
numInliers = maxNumInliers;

%fprintf('S1 time: %f\n', t1);
%fprintf('S2 time: %f\n', t2);
%fprintf('S3 time: %f\n', t3);

line.p1 = point1;
line.p2 = point2;
line.v = v;
line.inliers = inliers;
line.numInliers = numInliers;
line.percentInliers = numInliers / numPix;
end

function dists = pointLineDist(pos3d, v, p1)
    % http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    toPoint = repmat(p1, size(pos3d, 1), 1) - pos3d;
    v = repmat(v, size(pos3d, 1), 1);
    dists = cross(toPoint, v);
    dists = sqrt(sum(dists .^ 2, 2));
end
