classdef SuperPixUtils < Introspector
    properties(Constant)
        ID_PREFIX = 'id';
    end
    
    properties
        dh
        imset
        id
        
        superPix
        uniqueSupers
        maxSuperInd
        numSupers
        spixMap
        subsetMap
        loadedHash
        
        spixCache
        spixAvgPos
        spixAdj
        
        hline
        hlineFactor = 1;
        
        CACHE_FOLDER = 'spixCache';
    end
    
    methods
        function spixUtils = SuperPixUtils(dataHandler, imset, id)
            if nargin == 0
                return
            end
            self = spixUtils;
            self.dh = dataHandler;
            self.imset = imset;
            self.id = id;
            [self.superPix, self.uniqueSupers, self.numSupers] = ...
                        self.dh.getSpix(imset, id);
            cam = self.dh.getCam(self.imset, self.id);
            self.hline = cam.hline;
            self.maxSuperInd = max(self.uniqueSupers)+1;
            self.spixMap = uint16(zeros(self.maxSuperInd, 1));
            spixInds = 1:self.numSupers;
            self.spixMap(self.uniqueSupers+1) = spixInds;
        end
        
        function obj = constructDefault(~)
            obj = SuperPixUtils();
        end
        
        function cachePath = getCachePath(self, plusFile)
            fileId = '';
            if plusFile
                fileId = self.id;
            end
            props = self.dh.getStereoProps();
            basePath = fullfile(self.dh.getOutPath(self.imset), ...
                self.CACHE_FOLDER);
            cachePath = self.dh.genPath(basePath, 'cache', props, fileId);
        end
        
        function miny = getMinY(self)
            miny = floor(self.hline * self.hlineFactor);
        end
        
        function invalid = getInvalidSpix(self)
            spixBelowHline = self.superPix(self.getMinY():end, :);
            mask = ~ismember(self.uniqueSupers, spixBelowHline);
            invalid = self.uniqueSupers(mask);
        end
        
        function [uniqueSubset, numSubset, subsetHash] = getSpixSubset(self)
            spixBelowHline = self.superPix(self.getMinY():end, :);
            uniqueSubset = unique(sort(spixBelowHline));
            numSubset = numel(uniqueSubset);
            self.subsetMap = uint16(zeros(self.maxSuperInd, 1));
            spixInds = 1:numSubset;
            self.subsetMap(uniqueSubset+1) = spixInds;
            % DataHash function found at:
            % http://www.mathworks.com/matlabcentral/fileexchange/31272-datahash
            subsetHash = DataHash(uniqueSubset);
        end
        
        function result = isSubset(~, a, b)
            % Is "a" a subset of "b"
            result = all(ismember(a, b));
        end
        
        function inds = spixNumToIndex(self, nums, subset)
            if nargin < 3 || subset
                 % depends on getSpixSubset if subset is true
                if isempty(self.subsetMap)
                    self.getSpixSubset();
                end
                map = self.subsetMap;
            else
                map = self.spixMap;
            end
            inds = map(nums + 1);
        end
        
        function outStruct = iterSpix(self, in, args, func, opts)
            if nargin < 5
                opts = struct;
            end
            
            %allTic = tic;
            
            in.superPix = self.superPix;
            in.allUniqueSupers = self.uniqueSupers;
            in.allNumSupers = self.numSupers;
            if isfield(opts, 'spixNums')
                uniqueSubset = opts.spixNums;
                numSubset = numel(uniqueSubset);
                subsetHash = DataHash(uniqueSubset);
            else
                [uniqueSubset, numSubset, subsetHash] = self.getSpixSubset();
            end
            in.uniqueSupers = uniqueSubset;
            in.numSupers = numSubset;

            oneOnly = isfield(opts, 'one');
            if oneOnly;
                in.numSupers = 1;
                subsetHash = 0;
            end
            
            funcHasOut = false;
            if nargout(func)
                funcHasOut = true;
            end
            
            if isempty(self.spixCache) && ~oneOnly
                self.loadSpixCache();
            end
            cacheChanged = ~isequal(subsetHash, self.loadedHash);
            
            if isempty(self.spixAvgPos)
                self.spixAvgPos = zeros(self.maxSuperInd, 2);
            end
            
            function ind = sub2indFast(rows, row, col)
                % ~5x faster than sub2ind for 2d matrices
                ind = row + (col-1) * rows;
            end
            
            %funcTime = 0;
            %indTime = 0;
            %findTime = 0;
            %restTime = 0;
            spixIds = self.getSpixId(in.uniqueSupers);
            for i = 1:in.numSupers
                curSuper = in.uniqueSupers(i);
                spixId = spixIds{i};
                
                if isfield(self.spixCache, spixId)
                    spin = self.spixCache.(spixId);
                    y = spin.y; x = spin.x;
                    spin.pos2d = cat(2, x, y);
                    spin.spixIndex = sub2indFast(size(in.superPix, 1), y, x);
                    spin.spix = in.imShift(:, spin.spixIndex);
                else
                    %indTic = tic;
                    spixIndex2d = in.superPix == curSuper;
                    %indTime = indTime + toc(indTic);

                    %findTic = tic;
                    [y, x] = find(spixIndex2d);
                    % Linear indexing is MUCH faster than logical indexing so
                    % use that since we already use find
                    spixIndex = sub2indFast(size(spixIndex2d, 1), y, x);
                    avgPos2d = [mean(x(:)) mean(y(:))];
                    %findTime = findTime + toc(findTic);

                    %restTic = tic;
                    spin = struct;
                    spin.i = i;
                    spin.spixNum = curSuper;
                    spin.spixId = spixId;
                    spin.avgPos2d = avgPos2d;
                    spin.x = x;
                    spin.y = y;
                    spin.numPix = numel(x);
                    
                    self.spixCache.(spixId) = spin;
                    
                    % Vars that cannot be compressed in file easily but
                    % can be recreated easily so we don't store them
                    spin.spixIndex = spixIndex;
                    spin.pos2d = cat(2, x, y);
                    spin.spix = in.imShift(:, spixIndex);
                    %spin.spixIndex2d = spixIndex2d; % Takes up a lot of space
                    
                    %restTime = restTime + toc(restTic);
                    cacheChanged = true;
                end
                self.spixAvgPos(curSuper+1, :) = spin.avgPos2d;
                
                %functTic = tic;
                if funcHasOut
                    outStruct.(spixId) = func(in, args, spin);
                else
                    func(in, args, spin);
                end
                %funcTime = funcTime + toc(functTic);
            end
            if cacheChanged && ~oneOnly
                self.saveSpixCache(subsetHash);
            end
            %fprintf('Func time: %f\n', funcTime);
            %fprintf('Ind time: %f\n', indTime);
            %fprintf('Rest time: %f\n', restTime);
            %fprintf('Find time: %f\n', findTime);
            %fprintf('All time: %f\n', toc(allTic));
        end
        
        function saveSpixCache(self, subsetHash)
            cacheFolderPath = self.getCachePath(false);
            if ~exist(cacheFolderPath, 'dir')
                mkdir(cacheFolderPath);
            end
            cachePath = self.getCachePath(true);
            S = struct;
            S.spixCache = self.spixCache;
            S.subsetHash = subsetHash; %#ok<STRNU>
            % Update the hash representing what is loaded
            self.loadedHash = subsetHash; 
            save(cachePath, '-struct', 'S');
        end
        
        function savedHash = loadSpixCache(self)
            cachePath = self.getCachePath(true);
            if exist(cachePath, 'file')
                S = load(cachePath);
                self.spixCache = S.spixCache;
                savedHash = S.subsetHash;
                self.loadedHash = savedHash;
            else
                self.spixCache = struct;
                savedHash = 0;
            end
        end
        
        function emptyIterSpix(self)
            doNothing = @(x,y,z) 1;
            [~, in.imShift] = self.dh.getIm(self.imset, self.id);
            self.iterSpix(in, [], doNothing); 
        end
        
        function savedHash = loadAndUpdateCache(self, subsetHash)
            savedHash = 0;
            if isempty(self.spixCache)
                if nargin < 2
                    [~, ~, subsetHash] = self.getSpixSubset();
                end
                self.loadSpixCache();
            end
            if ~isequal(self.loadedHash, subsetHash)
                self.emptyIterSpix();
            end
        end
        
        function spixAvgPos = getSpixAvgPos(self)
            [~, ~, subsetHash] = self.getSpixSubset();
            self.loadAndUpdateCache(subsetHash);

            if isempty(self.spixAvgPos)
                %tic;
                spixAvgPos = zeros(max(self.uniqueSupers), 2);
                spixIds = fieldnames(self.spixCache);
                for i = 1:numel(spixIds)
                    spixId = spixIds{i};
                    spin = self.spixCache.(spixId);
                    spixAvgPos(spin.spixNum+1, :) = spin.avgPos2d;
                end
                self.spixAvgPos = spixAvgPos;
                %fprintf('Set avgPos: %f\n', toc);
            else
                spixAvgPos = self.spixAvgPos;
            end
        end
        
        function [spixAdj, borderTypes] = getSpixAdjacency(self)
            dispdir = self.dh.getStereoPath(self.imset);
            labelPath = fullfile(dispdir, sprintf('%s_label.txt', self.id));
            labelFile = fopen(labelPath, 'r');
            C = textscan(labelFile, '%d %d %d');
            lines = cat(2, C{:});
            [spixAdj, borderTypes] = getSpixAdjMatrix(lines);
            % Remove all invalid neighbours (above hline)
            invalid = self.getInvalidSpix() + 1;
            spixAdj(:, invalid) = -1;
            
            return;

            function [spixMat, borderTypes] = getSpixAdjMatrix(lines)
                maxSpixNum = max(lines(:)) + 1;
                spixMat = int16(zeros(maxSpixNum)) - 1;
                borderTypes = uint8(zeros(maxSpixNum)) - 1;
                for ti = 1:size(lines, 1)
                    line = lines(ti, :);
                    % The +1 is to avoid a 0 index which is invalid in
                    % matlab
                    spixNum = line(1) + 1;
                    spixNeighbour = line(2) + 1;
                    borderType = line(3);
                    spixMat(spixNum, spixNeighbour) = spixNeighbour;
                    borderTypes(spixNum, spixNeighbour) = borderType;
                    spixMat(spixNeighbour, spixNum) = spixNum;
                    borderTypes(spixNeighbour, spixNum) = borderType;
                end
            end
        end
        
        function neighbours = getNClosestNeighbours(self, spixNum, N)
            %tic;
            spixInd = spixNum + 1;
            if isempty(self.spixAdj)
                self.spixAdj = self.getSpixAdjacency();
            end
            around = self.spixAdj(spixInd, :);
            around = around(around > -1);
            level = around;
            while numel(around) < N + 1 ... % +1 because of self
                    && ~isempty(level) % shouldnt happen but prevents inf loop
                neighb = self.spixAdj(level, :);
                neighb = neighb(neighb > -1);
                level = unique(neighb(:));
                around = unique(cat(2, around, level'));
            end
            % Remove the target spix itself since it is a neighbour of its
            % immediate neighbours
            around = around(around ~= spixInd);
            
            avgPos = self.spixAvgPos;
            if isempty(avgPos)
                avgPos = self.getSpixAvgPos();
            end
            ourPos = avgPos(spixInd, :);
            aroundPos = avgPos(around, :);
            distPos = aroundPos - repmat(ourPos, size(aroundPos, 1), 1);
            %dists = sqrt(sum(distPos .^ 2, 2));
            dists = hypot(distPos(:, 1), distPos(:, 2));
            [~, ind] = sort(dists);
            sortedNeighbours = around(ind);
            neighbours = sortedNeighbours(1:N) - 1;
            %ntime = toc;
            %fprintf('Neighbour time: %f\n', ntime);
        end
    end
    
    methods(Static)
        function idPrefix = getIdPrefix()
            idPrefix = SuperPixUtils.ID_PREFIX;
        end
        
        function spixId = getSpixId(spixNums)
            if numel(spixNums) == 1
                % faster for single nums
                spixId = [SuperPixUtils.ID_PREFIX, sprintf('%d', spixNums)];
            else
                spixNumStrs = strtrim(cellstr(num2str(spixNums)));
                spixId = strcat(SuperPixUtils.ID_PREFIX, spixNumStrs);
            end
        end
        
        function spixNums = getNumsFromIds(spixIds)
            prefixEnd = length(SuperPixUtils.getIdPrefix()) + 1;
            spixNumStrs = char(spixIds);
            spixNumStrs = spixNumStrs(:, prefixEnd:end);
            spixNums = uint16(str2double(cellstr(spixNumStrs)));
        end
    end
end

