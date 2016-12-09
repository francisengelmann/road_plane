classdef RoadNet < GenericNet
    properties
        NET_SUBDIR = 'road';
        
        MIN_POS_PERCENT = 0.50;
        MIN_NEG_PERCENT = 0.99;
        
        numRegFeat;
        numContextFeat;
        contextSize = 4;
        useContextFeat = false;
        useContextWeights = false;
        
        % Postprocessing
        fitPlane = true;
        planeOpts = struct;
        bwConn = false;
        
        writeAngleMats = false;
        pitchMat;
        rollMat;
        planeOnly = false;
    end
    
    methods
        function roadNet = RoadNet(varargin)
            roadNet@GenericNet(varargin{:});
            roadNet.setTotalFeatures();
%           roadNet.planeOpts.maxDist = 0.15;
            roadNet.planeOpts.maxDist = 0.08;
            roadNet.planeOpts.maxNumPlanes = 1;
        end
        
        function obj = constructDefault(~)
            obj = RoadNet();
        end
        
        function totalFeat = setTotalFeatures(self, useContext)
            if nargin > 1
                self.useContextFeat = useContext;
            end
            [~, featVect] = self.getNumFeatures();
            numReg = numel(featVect);
            % Make the neighbour different (* 0.75) so weights arent NaN
            numCtxtFeat = numel(self.getContextFeatures(featVect .* 0.75, ...
                featVect));
            self.numRegFeat = numReg;
            self.numContextFeat = numCtxtFeat;
            totalFeat = numReg + numCtxtFeat;
            self.numFeatures = totalFeat;
            self.layers = totalFeat;
        end
        
        function exclude = getExcludeProps(self)
            exclude = {'pitchMat', 'rollMat'};
            exclude = cat(2, getExcludeProps@GenericNet(self), exclude);
        end
        
        function exclude = getExcludeChangedProps(self)
            exclude = {'writeAngleMats', 'planeOnly'};
            exclude = cat(2, getExcludeChangedProps@GenericNet(self), exclude);
        end
        
        function props = getFeatProps(self)
            props = {'useContextFeat', 'useContextWeights', ...
                'numContextFeat', 'contextSize'};
            props = cat(2, getFeatProps@GenericNet(self), props);
        end
        
        function props = getEvalProps(self)
            props = {'fitPlane', 'bwConn', 'planeOpts'};
            props = cat(2, getEvalProps@GenericNet(self), props);
        end
        
        function in = preprocessFeat(self, in)
            imsize = size(in.im);
            imsize = imsize(1:2);
            self.pitchMat = zeros(imsize) - 10;
            self.rollMat = zeros(imsize) - 10;
        end
        
        function [featureVect, featureStruct, in] = getFeatures(self, in, spin)
            boxSize = 10;
            
            im = in.im; pos = in.pos;
            spix = spin.spix; spixIndex = spin.spixIndex;

            imheight = size(im, 1);
            imwidth = size(im, 2);

            % The average x and y coordinates of the image are actually the most
            % important features it seems, sometimes simple is better
            avgx = spin.avgPos2d(1);
            avgy = spin.avgPos2d(2);

            % Horizon line
            hline = in.hline;
            belowHline = avgy > hline;

            %[boxPix, ~] = getImBox(im, avgx, avgy, boxSize);

            avgx = avgx / imwidth;
            avgy = avgy / imheight;
            % Distance from the middle of the image, the smaller the better, weighted
            % by y so that further down the image x can vary more than up.
            midDistW = (avgx - 0.5) / avgy;

            % 3D features
            if true
            pos3d = in.pos3d(:, spixIndex)';
            X = pos.X(spixIndex);
            Y = pos.Y(spixIndex);
            Z = pos.Z(spixIndex);

            % Try to fit a plane to the superpixel
            plane = ransacPlane(X, Y, Z);
            pnorm = plane.normal;
            % If the normal is pointing down flip it
            if pnorm(2) < 0
                pnorm = pnorm * -1;
            end
            % Get the % of pixels that are inliers to the plane
            percentInliers = plane.percentInliers;
            % Get the incline (angle relative to z axis) and twist (angle relative to
            % x axis) of the plane
            [pitch, roll] = getPlaneAngles(pnorm);
            if (max(spin.y) > 0.95 * imheight)
                % planes are pretty bad at the bottom of the image, just assume its
                % flat
                pitch = pitch * (pitch < 15);
                roll = roll * (roll < 15);
            end
            pitch = pitch / 90;
            roll = roll / 90;
            
            if self.writeAngleMats
                self.pitchMat(spixIndex) = pitch;
                self.rollMat(spixIndex) = roll;
            end

            % Median vs mean makes minimal difference
            avgPos = mean(pos3d);
            devPos = std(pos3d);
            end

            % Get the mean colour values, median here is worse
            meanRGB = mean(spix, 2) ./ 255;
            % Normalize the values so brightness does not affect them
            sumRGB = sum(meanRGB);
            normRGB = (meanRGB ./ sumRGB);
            devRGB = std(double(spix) ./ 255, false, 2);

            %featureVect = cat(2, meanRGB', avgx, avgy, normRGB');
            featureVect = cat(2, meanRGB', devRGB', normRGB', avgx, avgy, belowHline, ...
                pitch, roll, percentInliers, avgPos, devPos, ...
                midDistW);
                %midDistW, boxPix(:)');

            fs = struct;
            if true
            fs.meanRGB = meanRGB;
            fs.devRGB = devRGB;
            fs.avgx = avgx;
            fs.avgy = avgy;
            fs.belowHline = belowHline;
            fs.incline = pitch;
            fs.twist = roll;
            fs.percentInliers = percentInliers;
            fs.avgPos = avgPos;
            fs.devPos = devPos;
            fs.midDistW = midDistW;
            end
            featureStruct = fs;
        end
        
        function value = getFeatFromVect(~, vect, feats)
            if ~iscell(feats)
                feats = {feats};
            end
            inds = [];
            for i = 1:numel(feats)
                feat = feats{i};
                switch feat
                    case {'avgCol'}
                        ind = 1:3;
                    case {'devCol'}
                        ind = 4:6;
                    case {'pos2d'}
                        ind = 10:11;    
                    case {'pos3d'}
                        ind = 16:18;
                    case {'dev3d'}
                        ind = 19:21;
                    case {'angles'}
                        ind = 13:14;
                    otherwise
                        fprintf('unknown feature %s', feat);
                end
                inds = [inds ind]; %#ok<AGROW>
            end
            value = vect(:, inds);
        end
        
        function contextFeat = getContextFeatures(self, neighbFeat, refFeat)
            if ~self.useContextFeat
                contextFeat = [];
                return;
            end
            
            function res = wmean(vals, weights)
                weightsR = repmat(weights, 1, size(vals, 2));
                res = sum(vals .* weightsR, 1) ./ sum(weights);
            end
            
            featDiff = neighbFeat - repmat(refFeat, size(neighbFeat, 1), 1);
            
            diffPos2d = self.getFeatFromVect(featDiff, 'pos2d');
            % Make sure to include 1 as dim or mean will collapse all
            % features if there is only one neighbour
            avgDiffPos2d = mean(diffPos2d, 1);
            dists2d = hypot(diffPos2d(:,1), diffPos2d(:,2));
            avgDist = mean(dists2d);
            weights = 1 ./ (dists2d / min(dists2d));
            avgDiffFeats = {'avgCol', 'devCol', 'pos3d', 'dev3d', 'angles'};
            diffVals = self.getFeatFromVect(featDiff, avgDiffFeats);
            avgNeighbFeats = {'avgCol', 'devCol', 'pos2d', 'pos3d', ...
                'dev3d', 'angles'};
            neighbVals = self.getFeatFromVect(neighbFeat, avgNeighbFeats);
            if self.useContextWeights
                avgDiffAndNeighb = wmean([diffVals neighbVals], weights);
            else
                avgDiffAndNeighb = mean([diffVals neighbVals], 1);
            end
            stdNeighbFeats = {'avgCol', 'pos3d', 'angles'};
            stdNeighbVals = std(self.getFeatFromVect(neighbFeat, ...
                stdNeighbFeats), false, 1);
            contextFeat = [avgDiffAndNeighb avgDist avgDiffPos2d ...
                stdNeighbVals];
        end
        
        function features = getAllContextFeatures(self, features, spixUtils)
            if ~self.useContextFeat
                return;
            end
            
            tic;
            spixIds = fieldnames(features);
            spixNums = SuperPixUtils.getNumsFromIds(spixIds);
            featArr = cell2mat(struct2cell(features));
            spixMap = zeros(max(spixNums)+1, 1);
            spixInds = 1:numel(spixNums);
            spixMap(spixNums+1) = spixInds;
            
            %allFeat = zeros(numel(spixNums), self.numFeatures);
            for i = 1:numel(spixNums)
                spixId = spixIds{i};
                spixNum = spixNums(i);
                featInd = spixMap(spixNum+1);
                %featInd = spixUtils.spixNumToIndex(spixNum);
                featureVect = featArr(featInd, :);
                N = self.contextSize;
                neighbours = spixUtils.getNClosestNeighbours(spixNum, N);
                neighbInd = spixMap(neighbours+1);
                %neighbInd = spixUtils.spixNumToIndex(neighbours);
                neighbFeat = featArr(neighbInd, :);
                contextFeat = self.getContextFeatures(neighbFeat, featureVect);
                featWithContext = [featureVect contextFeat];
                %allFeat(i, :) = featWithContext;
                features.(spixId) = featWithContext;
                %features.(spixId) = contextFeat;
            end
            %features = cell2struct(num2cell(allFeat, 2), spixIds);
            fprintf('Context feature time: %f\n', toc);
        end
        
        function writePlaneAngleMatrices(self, imset, id)
            % Only write for testIds
            if self.writeAngleMats && any(ismember(self.testIds, id))
                outDir = self.getFeatPath(imset);
                pitchDir = fullfile(outDir, 'pitch');
                rollDir = fullfile(outDir, 'roll');
                if ~exist(pitchDir, 'dir')
                    mkdir(pitchDir);
                    mkdir(rollDir);
                end
                filename = [id '.txt'];
                sep = ' ';
                pitchPath = fullfile(pitchDir, filename);
                rollPath = fullfile(rollDir, filename);
                if self.OVERWRITE || ~exist(pitchPath, 'file')
                    dlmwrite(pitchPath, self.pitchMat, sep);
                    dlmwrite(rollPath, self.rollMat, sep);
                end
            end
        end
        
        function features = postprocessFeat(self, features, in, spixUtils)
            self.writePlaneAngleMatrices(in.imset, in.id);
            features = self.getAllContextFeatures(features, spixUtils);
        end
        
        function [imLabel, imConf, rest, ret] = postprocess(self, roadPix, ...
                roadConf, in, featureVects)
            spixAvgPos = self.getFeatFromVect(featureVects, 'pos3d');
            roadPixPos = spixAvgPos(in.spixLabel, :);
            allRoadSupers = in.labeledSupers;

            rest = struct;
            ret = struct;
            if self.fitPlane
                [roadPix, roadConf, mplane] = roadPlane(roadPix, roadConf, ...
                    roadPixPos, in.superPix, allRoadSupers, in.pos);
                rest.plane = mplane;
                rest.planes = planes;
                totalNumPlanes = numel(planes.planeCell);
                rest.numPlanes = totalNumPlanes;
                %out.allMinDists = allMinDists;
            end

            if self.bwConn
                CC = bwconncomp(roadPix);
                numPixels = cellfun(@numel, CC.PixelIdxList);
                [biggest, idx] = max(numPixels);
                roadPix(:) = false;
                roadPix(CC.PixelIdxList{idx}) = true;
            end

            imLabel = roadPix;
            if self.planeOnly
                imConf = false(size(roadConf));
            else
                imConf = roadConf;
            end
            
            function [roadPix, roadConf, mplane] = roadPlane(roadPix, roadConf, roadPixPos, superPix, allRoadSupers, pos)

                maxDistPlane = self.planeOpts.maxDist;

                planes = struct;
                %pks = [repmat(0.1, 1, 3), 0.02];
                pks = 0.02;
                %zsegs = {[-inf, 50], [50, 100], [100, inf], [-inf, inf]};
                zsegs = {[-inf, inf]};
                masks = {};
                maskStrs = {};
                for i = 1:numel(zsegs)
                    zseg = zsegs{i};
                    mask = (pos.Z >= zseg(1)) & (pos.Z < zseg(2));
                    masks = cat(1, masks, mask);
                    maskStrs = cat(1, maskStrs, ['z' num2str(abs(zseg(1))) ...
                        '_' num2str(zseg(2))]);
                end
                for i = 1:numel(masks)
                    pk = pks(i);
                    % Less chance to pick all 3 pixels within optimal plane
                    %pk = 0.02;
                    % Try and fit a plane to all the supposed road pixels
                    mask = masks{i} & roadPix;
                    plane = getPlane(pos, mask, maxDistPlane, pk);
                    planes.(maskStrs{i}) = plane;
                end
                mplane = plane;

                pk = 0.02;
                % Get a unique id for each road pixel so we can identify it after
                % using the inliers for linear indexing
                roadPixIds = reshape(1:numel(roadPix), size(roadPix));
                mask = roadPix;
                roadPixInliers = false(size(roadPix));
                MIN_INLIERS = numel(roadPix) * 0.0025;
                MAX_ROLL = 30;
                MAX_PITCH = 60;
                MIN_ZDEV = 0.5;
                MAX_LINE_INLIERPCENT = 0.70;
                planeCell = {};
                MAX_ITERS = 10;
                for i = 1:MAX_ITERS
                    if i > 1
                        plane = getPlane(pos, outliers2d, maxDistPlane, pk);
                    end
                    
                    try
                    [inliers2d, outliers2d] = ...
                        get2dInliers(roadPixIds, mask, plane.inliers);
                    catch
                        a=1;
                    end
                    
                    lineCell = {};
                    linePcentInliers = 0;
                    if i > 1
                        lineSumInliers = 0;
                        m2 = inliers2d;
                        MAX_LINES = 3;
                        for j = 1:MAX_LINES
                            if sum(m2(:)) < 3
                                break;
                            end
                            line = ransacLine(pos.X(m2), pos.Y(m2), ...
                                pos.Z(m2), maxDistPlane, pk);
                            lineCell = cat(1, lineCell, line);
                            [~, m2] = get2dInliers(roadPixIds, m2, ...
                                line.inliers);
                            lineSumInliers = lineSumInliers + line.numInliers;
                        end
                        linePcentInliers = lineSumInliers / plane.numInliers;
                    end
                    plane.lineCell = lineCell;

                    
                    [pitch, roll] = getPlaneAngles(plane.normal);
                    % zDev will be small for vertical planes
                    zDev = std(pos.Z(inliers2d));
                    if plane.numInliers > MIN_INLIERS && ...
                            abs(roll) < MAX_ROLL && abs(pitch) < MAX_PITCH ...
                            && zDev > MIN_ZDEV ...
                            && linePcentInliers < MAX_LINE_INLIERPCENT
                        plane.inliers2d = inliers2d;
                        planeCell = cat(1, planeCell, plane);
                        roadPixInliers = or(roadPixInliers, inliers2d);
                    end
                    
                    if sum(outliers2d(:)) < MIN_INLIERS ...
                            || numel(planeCell) >= self.planeOpts.maxNumPlanes
                        break;
                    end
                    
                    mask = outliers2d;
                end
                planes.planeCell = planeCell;
                if self.dh.DEBUG
                    self.showPlanes(planeCell, in);
                end

                % Also get all superpixels which have an average position close to the
                % plane
                toRoadPix = roadPixPos - ...
                    repmat(mplane.point, numel(allRoadSupers), 1);
                dists = abs(toRoadPix * mplane.normal');
                maxDistWeighted = maxDistPlane;
                % Let the height above the plane double in 50m
                %maxDistWeighted = maxDistPlane + roadPixPos(:, 3) / 50 * maxDistPlane;
                toKeep = dists <= maxDistWeighted;
                allRoadSupers = allRoadSupers(toKeep);
                roadPixSuper = ismember(superPix, allRoadSupers);

                %allMinDists = getMinPlaneDist(pos, planeCell);
                
                % Combine the inlier pixel and superpixel results
                roadPixPlane = or(roadPixInliers, roadPixSuper);
                removed = (roadPix - roadPixPlane) > 0;
                roadConf(removed) = roadConf(removed) * 0.5;
                roadPix = roadPixPlane;
            end
            
            function dists = getMinPlaneDist(pos, planeCell)
                allPos = [pos.X(:), pos.Y(:), pos.Z(:)];
                dists = inf(size(pos.X));
                for i = 1:numel(planeCell)
                    plane = planeCell{i};
                    toRoadPix = allPos - ...
                        repmat(plane.point, size(allPos, 1), 1);
                    thisDist = abs(toRoadPix * plane.normal');
                    thisDist = reshape(thisDist, size(pos.X));
                    dists = min(dists, thisDist);
                end
            end
            
            function dists = getMinPlaneDist2(pos, planeCell)
                allPos = [pos.X(:), pos.Y(:), pos.Z(:)];
                dists = inf(size(pos.X));
                for i = 1:numel(planeCell)
                    plane = planeCell{i};
                    if plane.normal(2) < 0
                        % make sure the normal is always facing up
                        plane.normal = plane.normal * -1;
                    end
                    toRoadPix = allPos - ...
                        repmat(plane.point, size(allPos, 1), 1);
                    thisDist = toRoadPix * plane.normal';
                    thisDist = reshape(thisDist, size(pos.X));
                    absLess = abs(thisDist) < abs(dists);
                    dists(absLess) = thisDist(absLess);
                end
            end
            
            function [inliers2d, outliers] = get2dInliers(ids2d, mask, inliers)
                % Get all the inlier pixels in terms of the road pixels
                maskedIds = ids2d .* mask;
                maksedLinear = maskedIds(maskedIds > 0);
                maskedInliers = maksedLinear(inliers);
                inliers2d = ismember(maskedIds, maskedInliers);
                outliers = mask & ~inliers2d;
            end
                        
            function plane = getPlane(pos, mask, varargin)
                plane = struct;
                if sum(mask(:)) < 3
                    yPos = -1.65;
                    if self.dh.COORD_TYPE == 1
                        yPos = yPos * -1;
                    end
                    plane.point = [0 yPos 0];
                    plane.u = [1 0 0];
                    plane.v = [0 0 1];
                    plane.normal = [0 1 0];
                    plane.inliers = false(size(mask));
                    plane.numInliers = 0;
                    plane.percentInliers = 0;
                    return;
                end

                plane = ransacPlane(pos.X(mask), pos.Y(mask), pos.Z(mask), ...
                        varargin{:});
                if plane.normal(2) < 0
                    plane.normal = plane.normal * -1;
                end
            end
        end
        
        function showPlanes(self, planeCell, in, exportDir, imPerPlane, ...
                plot2d, plotLines)
            im = in.im; cam = in.cam; pos = in.pos;
            if nargin < 4
                exportDir = false;
            end
            export = ~isequal(exportDir, false);
            if nargin < 5
                imPerPlane = true;
            end
            if nargin < 6
                plot2d = true;
            end
            if nargin < 7
                plotLines = false;
            end

            %X = pos.X; Y = pos.Y; Z = pos.Z;
            %postlims = {100, 100, 100};
            %postlims = {inf, inf, inf};
            lim = 150;
            [X, Y, Z] = self.dh.getXYZ(im, cam, pos.Z, lim);
            hline = cam.hline;
            cutoff = floor(hline * 0.7);
            ind = cutoff:size(im, 1);
            X = X(ind, :); Y = Y(ind, :); Z = Z(ind, :);

            numPlanes = numel(planeCell);
            cols = hsv(numPlanes+1);
            if imPerPlane
                outerIter = numPlanes;
            else
                outerIter = 1;
            end
            ALPHA = 0.5;
            for j = 1:outerIter
                f = getFig(1, export);
                set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
                imCol = double(im(ind, :, :))/255;
                %imCol = reshape(imCol, size(imCol, 1) * size(imCol, 2), 3);
                % Only plot every 5 pixels since matlab lags
                %ind = (1:5:numel(X));
                %h = scatter3(X(ind), Z(ind), Y(ind), 20, imCol(ind, :), 'filled');
                h = surf(X, Z, Y, 'cdata', imCol, 'EdgeColor', 'none');
                xlims = [min(X(:)) max(X(:))];
                zlims = [min(Z(:)) max(Z(:))];
                % Make the axis have the same scale
                axis equal;
                % The surface is hard to view unless you put some limits so it is centered
                xlim([-30 30]);
                % zlim is actually ylim
                viewYlim = [-5 20];
                zlim(viewYlim);
                if self.dh.COORD_TYPE == 1
                    set(gca,'ZDir','reverse');
                    zlim(fliplr(viewYlim * -1));
                end
                %alpha(h, 0.9);
                %view([0, 0]);
                hold on;

                if imPerPlane
                    innerStart = j;
                    innerIter = j;
                else
                    innerStart = 1;
                    innerIter = numPlanes;
                end
                ovIm = im;
                for i = innerStart:innerIter
                    plane = planeCell{i};
                    plotPlane(plane.point, plane.u, plane.v, ...
                        xlims, zlims, cols(i,:));
                    mask = plane.inliers2d;
                    scatter3(pos.X(mask), pos.Z(mask), pos.Y(mask), ...
                        10, cols(i+1,:));
                    if plot2d
                        ovIm = DataHandler.overlayImage(ovIm, mask, ...
                            cols(i+1,:), ALPHA);
                    end
                    if i > 1 && plotLines
                        for k = 1:numel(plane.lineCell)
                            line = plane.lineCell{k};
                            p1 = line.p1;
                            v = line.v;
                            mult1 = (zlims(2) - p1(3)) / v(3);
                            mult2 = (p1(3) - zlims(1)) / v(3);
                            linep = [p1 - v * mult2; p1; p1 + v * mult1];
                            plot3(linep(:,1), linep(:,3), linep(:,2), ...
                                'LineWidth',20);
                        end
                    end
                end
                if plot2d
                    f2 = getFig(2, export);
                    imshow(ovIm);
                end
                if export
                    jStr = num2str(j);
                    exportPath = fullfile(exportDir, [in.id '-p' jStr]);
                    export_fig(exportPath, f);
                    if plot2d
                        exPath2d = fullfile(exportDir, ...
                            [in.id '-2d' jStr]);
                        export_fig(exPath2d, f2);
                    end
                    
                end
            end
            hold off;
            
            function h = getFig(num, export)
                if export && ishandle(num)
                    set(0, 'CurrentFigure', num);
                    h = num;
                    clf(h);
                else
                    h = figure;
                end
            end
        end
        
        function multIds = getMultiPlaneIds(self, imset, ids)
            multIds = self.iterResults(imset, ids, @checkSingleId);
            function multiId = checkSingleId(in, res)
                if numel(res.planes.planeCell) > 1
                    multiId = in.id;
                else
                    multiId = '';
                end
            end
            emptyCells = cellfun(@isempty, multIds);
            % remove empty cells
            multIds(emptyCells) = [];
        end
        
        function iterShowPlanes(self, imset, ids, exportDir, varargin)
            if nargin < 4
                exportDir = 'resultsPlane';
            end
            exportParent = fullfile(self.dh.getOutSub(imset), exportDir);
            if ~exist(exportParent, 'dir')
                mkdir(exportParent);
            end
            opts.debug = true;
            self.iterResults(imset, ids, @handleSingleIm, opts);
            function out = handleSingleIm(in, res)
                out = [];
                f = figure(3);
                set(gcf, 'visible', 'off');
                ovIm = DataHandler.overlayImage(in.im, res.imLabel, ...
                            [1 0 0], 0.5);
                imshow(ovIm);
                exportPath = fullfile(exportParent, [in.id '-alabel']);
                export_fig(exportPath, f);
                close(f);
                self.showPlanes(res.planes.planeCell, in, exportParent, ...
                    varargin{:});
            end
        end
    end
    
    methods(Static)
        function dists = getPlaneDistLinear(posLin, plane, absv)
            toPlanePix = posLin - repmat(plane.point, size(posLin, 1), 1);
            dists = toPlanePix * plane.normal';
            if nargin < 3 || absv
                dists = abs(dists);
            end
        end
        
        function dists = getPlaneDist(pos, plane, varargin)
            allPos = [pos.X(:), pos.Y(:), pos.Z(:)];
            dists = RoadNet.getPlaneDistLinear(allPos, plane, varargin{:});
            dists = reshape(dists, size(pos.X));
        end
    end
end

