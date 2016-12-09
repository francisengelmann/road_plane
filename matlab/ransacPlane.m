function plane = ransacPlane(X, Y, Z, maxDist, pk, P)

if nargin < 6
    % Desired probability
    P = 0.99999;
end

if nargin < 5
    % Chance to pick all 3 pixels on the plane
    pk = 0.2;
end

if nargin < 4
    % Max 10 cm from plane
    maxDist = 0.1;
end

numPix = numel(X);
numSamples = round(log(1-P) / log(1-pk));
% 3 points for every sample
POINTS_PER_SAMP = 3;
replacement = numSamples * POINTS_PER_SAMP > numPix;
samples = randsample(numPix, numSamples * POINTS_PER_SAMP, replacement);
maxNumInliers = -1;
maxInliers = [];
bestAvgDist = inf;
% Return a xy plane by default
point = [0 0 0];
u = [1 0 0];
v = [0 1 0];
pnorm = [0 0 1];
% Get the 3d pos of each pixel in a numPixel x 3 array
% where the columns are (X,Y,Z)
pos3d = cat(2, X(:), Y(:), Z(:));
for i = 1:numSamples
    offset = (i-1) * POINTS_PER_SAMP;
    % Get the pixel number of each sample
    s1 = samples(offset+1);
    s2 = samples(offset+2);
    s3 = samples(offset+3);
    % Get the position of the pixel
    p1 = pos3d(s1, :);
    p2 = pos3d(s2, :);
    p3 = pos3d(s3, :);
    % Get the two vectors that represent the plane formed by the 3 points
    v1 = p2 - p1;
    v2 = p3 - p1;
    normal = cross(v1, v2);
    mag = norm(normal);
    if mag == 0
        % We chose parallel vectors, skip
        continue;
    end
    normal = normal ./ mag;
    
    % Since we normalized the normal the length of the projection is
    % just the dot product of the vector to the point and the normal
    % See: http://mathworld.wolfram.com/Point-PlaneDistance.html
    % and: http://mathworld.wolfram.com/Projection.html
    dists = pointPlaneDist(pos3d, p1, normal);
    inliers = dists <= maxDist;
    inDists = dists(inliers);
    avgDist = mean(inDists);
    % TODO Test squared distance instead as well
    %sqDist = sum(inDists .^ 2);
    numInliers = sum(inliers(:));
    if numInliers > maxNumInliers ...
            || (numInliers == maxNumInliers && avgDist < bestAvgDist)
        % If the inliers are the same and the distance is reduced, then
        % it is a better plane and use it
        maxNumInliers = numInliers;
        maxInliers = inliers;
        bestAvgDist = avgDist;
        point = p1;
        u = v1;
        v = v2;
        pnorm = normal;
    end
end
inliers = maxInliers;
numInliers = maxNumInliers;

plane.point = point;
plane.u = u;
plane.v = v;
plane.normal = pnorm;
plane.inliers = inliers;
plane.numInliers = numInliers;
plane.percentInliers = numInliers / numPix;

end

function dists = pointPlaneDist(pos3d, point, normal, sign)
    if nargin < 4
        sign = false;
    end
    toPoint = pos3d - repmat(point, size(pos3d, 1), 1);
    dists = toPoint * normal';
    if ~sign
        dists = abs(dists);
    end
end
