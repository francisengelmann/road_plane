classdef DataHandler < Introspector
    properties
        PROJECT_PATH
        SPSSTEREO_PATH
        STEREO_DIR
        DATA_DIR
        TRAIN_DIRS
        TEST_DIRS
        OUT_DIR
        OUT_SUB
        IM_LEFT = 'image_2';
        IM_RIGHT = 'image_3';
        IM_SUBDIRS = '0013';
        CALIB = 'calib';
        CALIB_FILE = '0013.txt';
        SET_TRAIN = 'road_training';  % training set of KITTI road benchmark
%       SET_TEST = 'stereo_multiview_training'; % KITTI object detection dataset
        SET_TEST = 'object_training';
        GT_TYPE = 'multiplane'
        numSpix = 1000;
        depthLims = struct;
        idDirMaps = struct;
        
        DEBUG = false;
        COORD_TYPE = 1;  % set to 1 so that Y is facing down
    end
    
    methods (Static, Access = protected)
        function getDefaults(dh)
            CODE_FOLDER = 'roadNet';
            % Go up until we reach the parent of the code folder
            dh.PROJECT_PATH = parentUntil(mfilename('fullpath'), CODE_FOLDER);
            dh.PROJECT_PATH = fullfile(dh.PROJECT_PATH, CODE_FOLDER);
            dh.SPSSTEREO_PATH = fullfile(dh.PROJECT_PATH, './spsstereo');
            dh.STEREO_DIR = 'results';
            dh.DATA_DIR = 'data/';
            dh.OUT_DIR = dh.DATA_DIR;
            dh.OUT_SUB = '';
            dh.depthLims.lim = 1000;
            dh.depthLims.postLims = {};
        end
    end
    
    methods
        function dh = DataHandler(outDir, dataDir, projectPath)
            DataHandler.getDefaults(dh);
            if nargin > 0
                dh.OUT_DIR = outDir;
            end
            if nargin > 1
                dh.DATA_DIR = dataDir;
            end
            if nargin > 2
                dh.PROJECT_PATH =  projectPath;
            end
            dh.DATA_DIR = pjoin(dh.PROJECT_PATH, dh.DATA_DIR);
            dh.OUT_DIR = pjoin(dh.PROJECT_PATH, dh.OUT_DIR);
            dh.SPSSTEREO_PATH = pjoin(dh.PROJECT_PATH, dh.SPSSTEREO_PATH);
            dh.TRAIN_DIRS = {pjoin(dh.DATA_DIR, dh.SET_TRAIN)};
            dh.TEST_DIRS = {pjoin(dh.DATA_DIR, dh.SET_TEST)};
        end
        
        function obj = constructDefault(~)
            obj = DataHandler();
        end
        
        function exclude = getExcludeProps(~)
            exclude = {'idDirMaps'};
        end
        
        function setDataDirs(self, trainDirs, testDirs)
            self.TRAIN_DIRS = getJoinedPaths(self, trainDirs);
            self.TEST_DIRS = getJoinedPaths(self, testDirs);
            
            function joinedDirs = getJoinedPaths(self, rawDirs)
                numDirs = numel(rawDirs);
                joinedDirs = cell(numDirs, 1);
                for i = 1:numDirs
                    newDir = pjoin(self.PROJECT_PATH, rawDirs{i});
                    joinedDirs{i} = newDir;
                end
            end
        end
        
        function dataDirs = getDataDirs(self, imset, id)
            if strcmp(imset, self.SET_TRAIN)
                dataDirs = self.TRAIN_DIRS;
            else
                dataDirs = self.TEST_DIRS;
            end
            if nargin > 2
                if ~isfield(self.idDirMaps, imset)
                    self.getIds(imset);
                end
                dirMap = self.idDirMaps.(imset);
                if ~isKey(dirMap, id)
                    dirNum = findIdDir(self, dataDirs, id);
                    dirMap(id) = dirNum;
                end
                dirNum = dirMap(id);
                dataDirs = dataDirs{dirNum};
            end
        end
        
        function dataDir = getDataDir(self, imset, id)
            dataDir = self.getDataDirs(imset, id);
        end
        
        function dirNum = findIdDir(self, dataDirs, id)
            found = false;
            for dirNum = 1:numel(dataDirs)
                dataDir = dataDirs{dirNum};
                imagePath = self.getImagePath(dataDir);
                image = strcat(imagePath, id, {'.png', '.jpg'});
                for i = 1:numel(image)
                    if exist(image{i}, 'file')
                        found = true;
                    end
                end
                if found
                    break;
                end
            end
            if ~found
                disp(id);
                error('Missing id');
            end
        end
        
        function imageDir = getImageDir(self, dataDir, right)
            if nargin < 3
                right = false;
            end
            
            subdir = self.IM_SUBDIRS;
            if iscell(subdir)
                % 1 if train 2 if test
                whichImset = 1 + any(strcmp(self.TEST_DIRS, dataDir));
                subdir = subdir{whichImset};
            end
            
            imDirs = {self.IM_LEFT, self.IM_RIGHT};
            % 1 for left, 2 for right
            index = right+1;
            imageDir = imDirs{index};
            if ~exist(fullfile(dataDir, imageDir), 'dir')
                imageList = dir(fullfile(dataDir, 'image*'));
                imageDir = imageList(index).name;
            end
            imageDir = fullfile(imageDir, subdir);
        end
        
        function [imagePath, imageDir] = getImagePath(self, dataDir, varargin)
            imageDir = self.getImageDir(dataDir, varargin{:});
            imagePath = fullfile(dataDir, imageDir);
        end
        
        function outPath = getOutPath(self, imset)
            outDir = self.OUT_DIR;
            whichImset = 1 + strcmp(imset, self.SET_TEST);
            if iscell(outDir)
                outPath = outDir{whichImset};
            else
                subdir = self.IM_SUBDIRS;
                if iscell(subdir)
                    subdir = subdir{whichImset};
                end
                outPath = fullfile(outDir, imset, subdir);
            end
        end
        
        function props = getStereoProps(~)
            props = {'numSpix'};
        end
        
        function stereoPath = getStereoPath(self, imset)
            props = self.getStereoProps();
            basePath = self.getOutPath(imset);
            stereoPath = self.genPath(basePath, self.STEREO_DIR, props);
        end
        
        function outSub = getOutSub(self, imset)
            outSub = fullfile(self.getOutPath(imset), self.OUT_SUB);
        end
        
        function imTransformPath = getImTransformPath(self, imset, id)
            if nargin > 2
                id = [id '.mat'];
            else
                id = '';
            end
            IM_TRANSF_DIR = 'image_transform';
            outPath = self.getOutPath(imset);
            imTransformPath = fullfile(outPath, IM_TRANSF_DIR, id);
        end
        
        function ids = getIds(self, imset)
            dataDirs = self.getDataDirs(imset);
            dirMap = containers.Map('KeyType', 'char', 'ValueType', 'uint8');
            ids = {};
            for dirNum = 1:numel(dataDirs)
                dataDir = dataDirs{dirNum};
                imagePath = self.getImagePath(dataDir);
                images = glob(imagePath, {'*.png', '*.jpg'});
                imageNames = {images.name};
                numFiles = size(imageNames, 2);
                fids = cell(numFiles, 1);
                for i = 1:numFiles
                    [~,name,~] = fileparts(imageNames{i});
                    fids{i} = name;
                end 
                ids = cat(1, ids, fids);
                dupKeyMask = isKey(dirMap, fids);
                if any(dupKeyMask)
                    dupKeys = fids(dupKeyMask);
                    disp(dupKeys);
                    error('Duplicate ids above');
                end
                curDirMap = containers.Map(fids, ...
                    repmat(uint8(dirNum), numel(fids), 1));
                dirMap = cat(1, dirMap, curDirMap);
            end
            self.idDirMaps.(imset) = dirMap;
        end
        
        function cam = getCam(self, imset, id)
            dataDir = self.getDataDir(imset, id);
            % read internal params
            calib_dir = fullfile(dataDir, self.CALIB);
            calibFile = self.CALIB_FILE;
            if iscell(calibFile)
                whichImset = 1 + strcmp(imset, self.SET_TEST);
                calibFile = calibFile{whichImset};
            end
            if ~isempty(calibFile)
                id = calibFile;
            end
            %calibPath = fullfile(calib_dir, sprintf('%s.txt', id(1:6)));
            calibPath = fullfile(calib_dir, '0019.txt');
            %disp(calibPath);
            %[K, ~, calib] = loadCalibration(calibPath);
            [~, K, calib] = loadCalibration(calibPath);
            
            %orig [~,~,tl] = KRt_from_P(calib.P_rect{3});  % left camera
            %orig [~,~,tr] = KRt_from_P(calib.P_rect{4});  % right camera
            [~,~,tl] = KRt_from_P(calib.P_rect{1});  % left camera
            [~,~,tr] = KRt_from_P(calib.P_rect{2});  % right camera
            f = K(1,1);
            baseline = abs(tr(1)-tl(1));   % distance between cams
            cam.f = f;
            cam.baseline = baseline;
            cam.K = K;
            % cam.P_left = calib.P_rect{3};
            %cam.P_right = calib.P_rect{4};
            cam.P_left = calib.P_rect{1};
            cam.P_right = calib.P_rect{2};
            cam.cu = K(1,3);
            cam.cv = K(2,3);
            cam.hline = K(2,3);
        end
        
        function cams = getCams(self, imset, ids)
            cams = struct;
            for i = 1:size(ids, 1)
                id = ids{i};
                sid = strcat('id', id);
                cams.(sid) = self.getCam(imset, id);
            end
        end
        
        function [im, imShift] = getIm(self, imset, id)
            dataDir = self.getDataDir(imset, id);
            imPath = self.getImagePath(dataDir);
            imfile = fullfile(imPath, sprintf('%s.png', id));
            if ~exist(imfile, 'file')
                imfile = fullfile(imPath, sprintf('%s.jpg', id));
            end
            im = imread(imfile);
            % Shift the image dimensions so we can index it easily
            imShift = shiftdim(im, 2);
        end
        
        function imTransforms = getImTransforms(self, imset, id, im)
            if nargin < 4
                im = self.getIm(imset, id);
            end
            
            imTransformPath = self.getImTransformPath(imset, id);
            function saveTforms(path)
                save(path, '-struct', 'imTransforms');
            end
            
            if exist(imTransformPath, 'file')
                imTransforms = load(imTransformPath);
            else
                imHsv = rgb2hsv(im);
                % Save space convert to uint8
                imHsv = uint8(imHsv .* 255);
                
                labCform = makecform('srgb2lab', ...
                    'AdaptedWhitePoint', whitepoint('D65'));
                %imLab = rgb2lab(im);
                imLab = applycform(im, labCform);
                
                lchCform = makecform('lab2lch');
                warning('off', 'images:encode_color:outputEncodingIgnored');
                imLch = applycform(imLab, lchCform);
                maxChroma = max(max(imLch(:, :, 2)));
                imLch = uint8(cat(3, imLch(:, :, 1) / 100, ...
                    imLch(:, :, 2) / maxChroma, imLch(:, :, 3) / 360) * 255);
                
                imTransforms.hsv = imHsv;
                imTransforms.lab = imLab;
                imTransforms.lch = imLch;
                try
                    saveTforms(imTransformPath);
                catch ME
                    if strcmp(ME.identifier, 'MATLAB:save:noParentDir')
                        parentDir = fileparts(imTransformPath);
                        mkdir(parentDir);
                        saveTforms(imTransformPath);
                    else
                        rethrow(ME);
                    end
                end
            end
        end
        
        function [pos, neg, none, gtRaw] = getGroundTruth(self, id)
            imset = self.SET_TRAIN;
            dataDir = self.getDataDir(imset, id);
            if isempty(strfind(id, '_'))
                gtRoadName = ['road_', id];
            else
                gtRoadName = strrep(id, '_', '_road_');
            end
            gtFilename = [gtRoadName '.png'];
            imDir = self.getImageDir(dataDir);
            gtImDir = ['gt_', imDir];
            gtRoadPath = fullfile(dataDir, gtImDir, gtFilename);
            gt = imread(gtRoadPath);
            if strcmp(self.GT_TYPE, 'binary') && size(gt, 3) > 1
                % If forcing binary then just take the positive plane from
                % the new gt images
                gt = uint8(gt(:, :, 1) & gt(:, :, 3)) * 255;
            end
            if size(gt, 3) > 1
                % New multi-plane ground truth
                pos = gt(:, :, 1) & gt(:, :, 3);
                neg = gt(:, :, 1) & ~gt(:, :, 3);
                none = ~(pos | neg); % or: ~(any(gt, 3))
            else
                % Simple binary (black/white) ground truth
                pos = gt > 0;
                neg = ~pos;
                none = false(size(gt));
            end
            gtRaw = gt;
        end
        
        function [spix, uniqueSupers, numSupers] = getSpix(self, imset, id)
            dispdir = self.getStereoPath(imset);
            spfile = fullfile(dispdir, sprintf('%s_segment.png', id)); 
            if ~exist(spfile, 'file')
                error('you haven''t ran spsstereo code yet for: %s...\n', id);
            end
            spix = imread(spfile);
            spix = uint16(spix);
            uniqueSupers = unique(sort(spix));
            numSupers = numel(uniqueSupers);
        end
        
        function disp = getDisparity(self, imset, id)
            origWD = pwd;
            outdir = self.getStereoPath(imset);
            if ~exist(outdir, 'dir')
                mkdir(outdir);
            end;
            outfile = fullfile(outdir, sprintf('%s_left_disparity.png', id)); 
            dataDir = self.getDataDir(imset, id);
            imPathLeft = self.getImagePath(dataDir);
            imPathRight = self.getImagePath(dataDir, true);
            imFilename = sprintf('%s.png', id);
            imfile_left = fullfile(imPathLeft, imFilename);
            imfile_right = fullfile(imPathRight, imFilename);
            spsstereo = fullfile(self.SPSSTEREO_PATH, 'spsstereo');
            if ~exist(outfile, 'file')
                cd(outdir);
                c = onCleanup(@() cd(origWD));
                fprintf('running spsstereo, may take a few secs...\n');
                tic;
                clean = false;
                if ~exist(imfile_left, 'file')
                    clean = true;
                    convertToPng(imfile_left)
                    convertToPng(imfile_right)
                end
                cmd = sprintf('"%s" "%s" "%s" %d', spsstereo, ...
                    imfile_left, imfile_right, self.numSpix);
                system(cmd);
                e=toc;
                if clean
                    cleanup(imfile_left)
                    cleanup(imfile_right)
                end
                fprintf('finished! total time: %0.4f\n', e);
                if ispc()
                    % If we are on PC it doesn't output in the working
                    % directory for some reason, so move it
                    fileWildcard = sprintf('%s_*', id);
                    source = fullfile(imPathLeft, fileWildcard);
                    movefile(source, outdir);
                end
            end;
            sprintf('Saving disparity: %s',outfile);

            disp = imread(outfile);
            disp = double(disp)/256;

            function convertToPng(imfile)
                imfile_or = strrep(imfile, '.png', '.jpg');
                im = imread(imfile_or);
                imwrite(im, imfile);
            end

            function cleanup(imfile)
                delete(imfile);
            end
        end
        
        function Z = getDepth(self, imset, id, cam)
            disp = self.getDisparity(imset, id);
            if self.COORD_TYPE == 1
                camdata = cam;
                disparity = disp;
                disparity = single(disparity);
                disparity(disparity == 0) = 0.1;
                Z = camdata.f * camdata.baseline ./ double(disparity);
            else
                f = cam.f;
                Z = f * cam.baseline ./ disp;
                lim = 1e6;
                Z = min(Z, lim);
            end
        end
        
        function [X, Y, Z] = getXYZ(self, im, cam, depth, limit, postLims)
            if nargin < 5
                limit = inf;
            end
            if nargin < 6 || isempty(postLims)
                postLims = {inf, inf, inf};
            end

            Z = depth;
            Z = min(Z, limit);
            if self.COORD_TYPE == 1
                camdata = cam;
                RD = Z;
                [xx, yy] = meshgrid(1:size(RD,2), 1:size(RD,1));
                u = (xx(:) - camdata.cu) .* RD(:) / camdata.f;
                v = (yy(:) - camdata.cv) .* RD(:) / camdata.f;
                xyz = [u(:), v(:), RD(:)];
                X = reshape(xyz(:, 1), size(Z));
                Y = reshape(xyz(:, 2), size(Z));
            else
                x = size(im, 2);
                y = size(im, 1);

                f = cam.f;
                px_d = cam.K(1,3);
                py_d = cam.K(2,3);
                [xvals, yvals] = meshgrid(1:x, y:-1:1);
                % Calculate X and Y from depth
                X = (xvals - px_d) .* (Z / f);
                Y = (yvals - py_d) .* (Z / f);
            end                

            X = applyLims(X, postLims{1});
            Y = applyLims(Y, postLims{2});
            Z = applyLims(Z, postLims{3});
            
            function var = applyLims(var, lims)
                if max(size(lims)) > 1
                    var = min(lims(2), max(var, lims(1)));
                else
                    var = min(lims, max(var, -lims));
                end
            end
        end

        function [outStruct, outCell, totalTime] = iterIds(self, imset, ids, opts, func, args)
            if nargin < 6
                args = 0;
            end

            numIds = numel(ids);
            sids = strcat('id', ids);

            funcHasOut = false;
            if nargout(func)
                funcHasOut = true;
            end
            
            if ~isfield(opts, 'cellOnly')
                emptyStruct = struct;
                if funcHasOut
                    emptyStruct = func(0, 0, true);
                end 
                outStructArr(numIds) = emptyStruct;
                cellOnly = false;
            else
                outCell = cell(numIds, 1);
                cellOnly = true;
            end

            totalTic = tic;
            
            if self.DEBUG == true || isfield(opts, 'debug')
                M = 0;
            else
                M = inf;
            end
            %parfor (i = 1:numIds, M)
            for (i = 1:numIds)
                iterTic = tic;

                id = ids{i};
                %underscorePos = strfind(id, '_');
                %if isempty(underscorePos)
                %    idNum = id;
                %else
                %    idNum = id(underscorePos(end)+1:end);
                %end
                %idNum = str2double(idNum);
                %rng(idNum*10);

                in = struct;
                in.i = i;
                in.imset = imset;
                
                disp('DataHandler.m line 542 FrameId: ');
                disp(id);
                in.id = id;
                sid = sids{i};
                in.sid = sid;
                
                if isfield(opts, 'prePrint')
                    fprintf('Starting: %s;  i: %d\n', id, i);
                end
                
                if ~isfield(opts, 'noIm')
                    [in.im, in.imShift] = self.getIm(imset, id); %#ok<PFBNS>
                    %in.imTransforms = self.getImTransforms(imset, id, in.im);
                end

                if isfield(opts, 'groundTruth')
                    [in.gtPos, in.gtNeg, in.gtNone, in.gtRaw] = ...
                        self.getGroundTruth(id);
                end
                
                % this is likely where the computation happens
                if ~isfield(opts, 'noPos')
                    cam = self.getCam(imset, id);
                    in.hline = cam.hline;
                    Z = self.getDepth(imset, id, cam);
                    in.depth = Z;
                    if isfield(opts, 'depthLims')
                        lim = opts.depthLim;
                        postLims = opts.postLims;
                    else
                        lim = self.depthLims.lim;
                        postLims = self.depthLims.postLims;
                    end
                    [X, Y, Z] = self.getXYZ(in.im, cam, Z, lim, postLims);
                    pos = struct;
                    pos.X = X; pos.Y = Y; pos.Z = Z;
                    pos3d = shiftdim(cat(3, X, Y, Z), 2);
                    in.cam = cam;
                    in.pos = pos;
                    in.pos3d = pos3d;
                end

                if isfield(opts, 'spix')
                    [in.superPix, in.uniqueSupers, in.numSupers] = ...
                        self.getSpix(imset, id);
                end
                
                % this is where it crashes
                if funcHasOut
                    funcOut = func(in, args, false);  % #ok<PFBNS>
                else
                    func(in, args, false);
                    funcOut = struct;
                end

                if ~cellOnly
                    outStructArr(i) = funcOut;
                else
                    outCell{i} = funcOut;
                end

                iterTime = toc(iterTic);
                if ~isfield(opts, 'noPrint')
                    fprintf('Processing time for: %s = %f\n', id, iterTime);
                end
            end

            if ~cellOnly
                [outStruct, outCell] = ...
                    DataHandler.reorgStruct(outStructArr, sids);
            else
                outStruct = outCell;
            end

            totalTime = toc(totalTic);
            if ~isfield(opts, 'noPrint')
                fprintf('Total time: %f\n', totalTime);
            end
        end

    end
    
    methods(Static)
        function [trainIds, testIds, randOrd] = splitIds(ids)
            numIds = numel(ids);
            randOrd = randperm(numIds);

            half = floor(numIds/2);
            trainIds = ids(randOrd(1:half));
            testIds = ids(randOrd(half+1:end));

            trainIds = sort(trainIds);
            testIds = sort(testIds);
        end
        
        function [outStruct, outCell] = reorgStruct(outStructArr, sids)
            outCell = squeeze(struct2cell(outStructArr));
            if size(outCell, 2) == 1
                outCell = outCell';
            end
            fields = fieldnames(outStructArr);
            outStruct = struct;
            for j = 1:numel(fields)
                field = fields{j};
                outStruct.(field) = cell2struct(outCell(j, :)', sids);
            end
        end
        
        function ovIm = overlayImage(im, mask, col, alpha)
            im = im2double(im);
            maskRGB = repmat(mask, [1 1 3]);
            untocuhed = double(maskRGB == false) .* im;
            col = alpha * col;
            colArr = repmat(reshape(col, 1, 1, 3), [size(mask) 1]);
            overlayComponent = maskRGB .* colArr;
            origImageComponent = maskRGB .* (1 - alpha) .* im;
            ovIm = untocuhed + overlayComponent + origImageComponent;
        end
        
        function paths = pjoin(paths1, paths2)
            paths = {};
            numelPaths1 = max(1, numel(paths1) * iscell(paths1));
            numelPaths2 = max(1, numel(paths2) * iscell(paths2));
            for i = 1:numelPaths1
                if iscell(paths1)
                    path1 = paths1{i};
                else
                    path1 = paths1;
                end
                for j = 1:numelPaths2
                    if iscell(paths2)
                        path2 = paths2{j};
                    else
                        path2 = paths2;
                    end
                    if path2(1) == '/' || (numel(path2) > 1 && path2(2) == ':')
                        % Path 2 is absolute, do not join paths
                        path = path2;
                    else
                        path = fullfile(path1, path2);
                    end
                    paths = cat(1, paths, path);
                end
            end
            if numel(paths) == 1
                % Single path, don't return a cell
                paths = paths{1};
            end
        end
    end
end

function paths = pjoin(paths1, paths2)
    paths = DataHandler.pjoin(paths1, paths2);
end

function fileNames = glob(path, wildcards)
    for i = 1:numel(wildcards)
        newNames = dir(fullfile(path, wildcards{i}));
        if i == 1
            fileNames = newNames;
        else
            fileNames = cat(1, fileNames, newNames);
        end
    end
end

function [path, pathWith] = parentUntil(startPath, until)
    [path, name, ~] = fileparts(startPath);
    while ~strcmp(name, until)
        [path, name, ~] = fileparts(path);
    end
    pathWith = fullfile(path, name);
end
