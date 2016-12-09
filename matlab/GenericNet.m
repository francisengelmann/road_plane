classdef (Abstract) GenericNet < Introspector
    %#ok<*PROP>
    properties
        % DataHandler
        dh;
        
        inputs;
        targets;
        net;
        tr;
        accuracy;
        
        trainIds;
        testIds;
        
        NET_DIR = 'roadnet';
        FEAT_DIR = 'feat';
        TRAIN_DATA_DIR = 'trained';
        TRAIN_DATA_FNAME = 'trainData';
        TRAINED_NET_DIR = 'trainedNet';
        TRAINED_NET_FNAME = 'trainedNet';
        RESULT_DIR = 'eval';
        IMLABEL_DIR = 'resultsSeg';
        
        numFeatures;
        featVect;
        layers;
        epochs = 40000;
        maxFail = 100;
        ratios = [.75 .20 .05];
        useGPU = 'only';
        trainFcn = 'trainscg';
        mapstd = false;
        
        THRESH = 0.5;
        
        OVERWRITE = false;
        NESTED_OUT = true;
        EVAL_TEST_PROP = '';
        
        existNetPath = '';
        testOnly = false;
    end
    
    properties (Abstract)
        NET_SUBDIR;
        
        % Minimum positive percent of ground truth for it to be considered
        % a positive sample
        MIN_POS_PERCENT;
        % Similarly minium percent for a sample to be considered negative
        MIN_NEG_PERCENT;
    end
    
    methods (Abstract)
       [featureVect, featureStruct] = getFeatures(self, in, spin);
    end
    
    methods
        function nnet = GenericNet(dataHandler, trainIds, testIds, existingNetPath)
            if nargin < 1 || isequal(dataHandler, false)
                dataHandler = DataHandler();
            end
            nnet.dh = dataHandler;
            
            if nargin < 4
                existingNetPath = '';
            end
            if ~isempty(existingNetPath)
                existingNetPath = DataHandler.pjoin(...
                    dataHandler.PROJECT_PATH, existingNetPath);
                nnet.testOnly = true;
                nnet.disableInspect = true;
            end
            nnet.existNetPath = existingNetPath;
            
            if (nargin < 2 || ~iscell(trainIds)) && ~nnet.testOnly
                trainIds = dataHandler.getIds(dataHandler.SET_TRAIN);
            end
            if nargin < 3 || ~iscell(testIds)
                testIds = dataHandler.getIds(dataHandler.SET_TEST);
            end
            nnet.trainIds = trainIds;
            nnet.testIds = testIds;
        end
        
        function [featVect, numFeat] = initNumFeatures(self, setLayers)
            if nargin < 2
                setLayers = false;
            end
            if ~self.testOnly
                featVect = self.getSingleFeatVect();
            else
                % If created just for eval then get num feats from loading net
                self.loadNet();
                featVect = zeros(1, self.net.input.size);
            end
            self.featVect = featVect;
            numFeat = numel(featVect);
            self.numFeatures = numFeat;
            if setLayers
                self.layers = self.numFeatures;
            end
        end
        
        function [numFeat, featVect] = getNumFeatures(self)
            if isempty(self.numFeatures) || isempty(self.featVect)
                self.initNumFeatures();
            end
            numFeat = self.numFeatures;
            featVect = self.featVect;
        end
        
        function exclude = getExcludeProps(~)
            exclude = {'inputs', 'targets', 'net', 'tr'};
        end
        
        function exclude = getExcludeChangedProps(~)
            exclude = {'trainIds', 'testIds', 'accuracy', 'EVAL_TEST_PROP', ...
                'featVect'};
        end
        
        function include = getIncludeProps(~)
            include = {'dh.numSpix', 'dh.depthLims.lim', ...
                'dh.depthLims.postLims'};
        end
        
        function netSubdir = getNetSubdir(self, imset)
            netSubdir = self.NET_SUBDIR;
            if iscell(netSubdir)
                whichImset = 1 + strcmp(imset, self.SET_TEST);
                netSubdir = netSubdir{whichImset};
            end
        end
        
        function netPath = getNetPath(self, imset)
            netPath = fullfile(self.dh.getOutPath(imset), self.NET_DIR, ...
                self.getNetSubdir(imset));
        end
        
        function props = getFeatProps(~)
            props = {'numFeatures', 'dh.numSpix', 'dh.depthLims.lim', ...
                'dh.depthLims.postLims', 'dh.COORD_TYPE'};
        end
        
        function featPath = getFeatPath(self, imset, id)
            if nargin < 3
                id = '';
            end
            
            netPath = self.getNetPath(imset);
            if ~self.testOnly
                self.getNumFeatures();
                props = self.getFeatProps();
            else
                props = {};
            end
            featPath = self.genPath(netPath, self.FEAT_DIR, props, id);
        end
        
        function props = getTrainDataProps(~)
            props = {'MIN_NEG_PERCENT', 'MIN_POS_PERCENT', 'dh.GT_TYPE'};
        end
        
        function trainDataPath = getTrainDataPath(self, imset, plusFile)
            if nargin < 2 || isempty(imset);
                imset = self.dh.SET_TRAIN;
            end
            if nargin < 3 || ~plusFile
                fileName = '';
            else
                fileName = self.TRAIN_DATA_FNAME;
            end
            props = self.getTrainDataProps();
            m100 = @(x) x .* 100;
            ops = {m100, m100, []};
            featPath = self.getFeatPath(imset);
            trainDataPath = self.genPath(featPath, self.TRAIN_DATA_DIR, ...
                props, fileName, ops);
        end
        
        function props = getTrainedNetProps(~)
            props = {'layers', 'maxFail', 'ratios', 'trainFcn'};
        end
        
        function trainedNetPath = getTrainedNetPath(self, imset, plusFile)
            if nargin < 2 || isempty(imset)
                imset = self.dh.SET_TRAIN;
            end
            if nargin < 3 || ~plusFile
                fileName = '';
            else
                fileName = self.TRAINED_NET_FNAME;
            end
            props = self.getTrainedNetProps();
            ops = cell(4, 1);
            ops{3} = @(x) x .* 100;
            trainDataPath = self.getTrainDataPath(imset);
            trainedNetPath = self.genPath(trainDataPath, ...
                self.TRAINED_NET_DIR, props, fileName, ops);
        end
        
        function props = getEvalProps(~)
            props = {};
        end
        
        function resultsPath = getResultsPath(self, imset, id)
            if nargin < 3
                id = '';
            end
            if 0 %~self.testOnly
                props = self.getEvalProps();
                tnetPath = self.getTrainedNetPath(imset);
            else
                % testing only, produce a non-nested path
                props = {};
                tnetPath = self.getNetPath(imset);
            end
            resultsPath = self.genPath(tnetPath, self.RESULT_DIR, props, id);
        end
        
        function outPath = getOutSubPath(self, imset, subdir, labelsOnly)
            if nargin < 3 || isempty(subdir)
                subdir = self.IMLABEL_DIR;
            end
            if nargin < 4
                labelsOnly = false;
            end
            if ischar(labelsOnly)
                subdir = [subdir labelsOnly];
            elseif labelsOnly
                subdir = [subdir 'Bin'];
            end
            if self.NESTED_OUT
                resultsPath = self.getResultsPath(imset);
                outPath = resultsPath;
            else
                outPath = self.dh.getOutSub(imset);
            end
            outPath = fullfile(outPath, subdir);
        end
        
        function swapIdSets(self)
            tmp = self.trainIds;
            self.trainIds = self.testIds;
            self.testIds = tmp;
        end
        
        function featureVect = getSingleFeatVect(self)
            if ~isempty(self.trainIds)
                imset = self.dh.SET_TRAIN;
                id = self.trainIds(1);
            else
                imset = self.dh.SET_TEST;
                id = self.testIds(1);
            end
            sid = strcat('id', id{1});
            
            opts.noPrint = true;
            opts.debug = true;
            outStruct = self.dh.iterIds(imset, id, opts, @getVect);
            function out = getVect(in, args, getStruct)
                if getStruct
                    out.featureVect = 0;
                    return;
                end
                in = self.preprocessFeat(in);
                sopts.one = true;
                features = [];
                spixUtils = SuperPixUtils(self.dh, in.imset, in.id);
                spixUtils.iterSpix(in, args, @handleSpix, sopts);
                function t = handleSpix(in, ~, spin)
                    t = 0;
                    features = self.getFeatures(in, spin);
                end
                %features = self.postprocessFeat(features);
                out.featureVect = features;
            end
            featureVect = outStruct.featureVect.(sid);
        end
        
        function in = preprocessFeat(~, in)
        end
        
        function prepFeatures(self, imset, ids, overwrite)
            if nargin < 4
                overwrite = false;
            end
            
            featPath = self.getFeatPath(imset);
            if ~exist(featPath, 'dir')
                mkdir(featPath);
            end
            
            opts = struct;
            self.dh.iterIds(imset, ids, opts, @prepSingleIm);
            
            function prepSingleIm(in, args, ~)                   
                if ~overwrite && ~self.OVERWRITE && ...
                        exist(self.getFeatPath(imset, in.id), 'file')
                    return;
                end
                
                in = self.preprocessFeat(in);
                
                spixUtils = SuperPixUtils(self.dh, in.imset, in.id);
                features = spixUtils.iterSpix(in, args, @handleSpix);
                
                function featureVect = handleSpix(in, ~, spin)
                    featureVect = self.getFeatures(in, spin);
                end
                
                features = self.postprocessFeat(features, in, spixUtils);
                
                self.saveFeatures(imset, in.id, features)
            end
        end
        
        function features = postprocessFeat(~, features, in, spixUtils)
        end
        
        function saveFeatures(self, imset, id, features)
            featPath = self.getFeatPath(imset, id);
            save(featPath, 'features');
        end
        
        function features = loadFeatures(self, imset, id)
            featPath = self.getFeatPath(imset, id);
            S = load(featPath);
            features = S.features;
        end
        
        function [inputs, targets] = prepTrainData(self, overwrite)
            if nargin < 2
                overwrite = false;
            end
            
            if ~overwrite && ~self.OVERWRITE ...
                    && exist(self.getTrainDataPath('', true), 'file')
                return;
            end
            
            imset = self.dh.SET_TRAIN;
            ids = self.trainIds;
            
            tic;

            opts.groundTruth = true;
            outStruct = self.dh.iterIds(imset, ids, opts, @prepSingleIm);
            
            function out = prepSingleIm(in, args, getStruct)
                if nargin > 2 && getStruct
                    out.samples = 0; out.targets = 0;
                    return;
                end

                % Only one class, road or not
                numClasses = 1;

                spixUtils = SuperPixUtils(self.dh, in.imset, in.id);
                samples = zeros(spixUtils.numSupers, self.getNumFeatures());
                truth = zeros(numClasses, spixUtils.numSupers);
                numSamples = 0;
                
                features = self.loadFeatures(imset, in.id);

                spixUtils.iterSpix(in, args, @handleSpix);
                
                function sout = handleSpix(in, ~, spin)    
                    sout = 0;
                    
                    numPix = spin.numPix;
                    gtPos = in.gtPos(spin.spixIndex);
                    gtNeg = in.gtNeg(spin.spixIndex);
                    percentPos = sum(gtPos(:)) / numPix;
                    percentNeg = sum(gtNeg(:)) / numPix;

                    positiveExample = percentPos >= self.MIN_POS_PERCENT;
                    negativeExample = percentNeg >= self.MIN_NEG_PERCENT;
                    if positiveExample || negativeExample            
                        numSamples = numSamples + 1;

                        featureVect = features.(spin.spixId);
                        samples(numSamples, :) = featureVect;
                        truth(numSamples) = positiveExample;
                    end
                end

                % prune any preallocated that were not used
                samples = samples(1:numSamples, :);
                truth = truth(1:numSamples);

                out.samples = samples;
                out.targets = truth;
            end
            
            inputs = cell2mat(struct2cell(outStruct.samples));
            targets = cell2mat(struct2cell(outStruct.targets)');
            self.inputs = inputs;
            self.targets = targets;
            
            self.saveTrainData(inputs, targets);

            prepTime = toc;
            fprintf('Prep time: %f\n', prepTime);
        end
        
        function saveTrainData(self, inputs, targets)
            trainDataPath = self.getTrainDataPath('', true);
            trainDataDir = fileparts(trainDataPath);
            if ~exist(trainDataDir, 'dir')
                mkdir(trainDataDir);
            end
            save(trainDataPath, 'inputs', 'targets');
        end
        
        function loadTrainData(self)
            trainDataFile = self.getTrainDataPath('', true);
            S = load(trainDataFile);
            self.inputs = S.inputs;
            self.targets = S.targets;
        end
        
        function train(self, overwrite)
            if nargin < 2
                overwrite = false;
            end 
            
            if ~overwrite && ~self.OVERWRITE ...
                    && exist(self.getTrainedNetPath('', true), 'file')
                return;
            end
            
            %rng(1);
            
            if isempty(self.inputs)
                self.loadTrainData();
            end

            % pattern recognition network
            hiddenLayerSize = self.layers;
            net = patternnet(hiddenLayerSize);
            net.trainFcn = self.trainFcn;
            if self.mapstd
                net.inputs{1}.processFcns = [net.inputs{1}.processFcns, 'mapstd'];
            end
            net.trainParam.epochs = self.epochs;
            net.trainParam.max_fail = self.maxFail;

            % divide data for training, validation and testing
            net.divideParam.trainRatio = self.ratios(1);
            net.divideParam.valRatio = self.ratios(2);
            net.divideParam.testRatio = self.ratios(3);

            % train
            tic;
            useParallel = 'no';
            if gpuDeviceCount() == 0
                useParallel = 'yes';
            end
            [net, tr] = train(net, self.inputs', self.targets, ...
                'useGPU', self.useGPU, 'useParallel',useParallel);
            trainTime = toc;
            fprintf('Train time: %f\n', trainTime);
            outputs = net(self.inputs');
            errors = abs(gsubtract(self.targets, outputs));
            performance = perform(net, self.targets, outputs);
            self.net = net;
            self.tr = tr;

            % visualize
            %figure(1), plotperform(tr);
            %figure(2), plottrainstate(tr);
            %figure(3), plotconfusion(targets, outputs);
            %figure(4), ploterrhist(errors);

            % a simple way to compute accuracy
            diff = (outputs > 0.5) - self.targets;
            self.accuracy = sum(diff == 0)/size(vec2ind(self.targets)',1);

            self.saveTrainedNet(net, tr, self.accuracy);
        end
        
        function saveTrainedNet(self, net, tr, accuracy)
            trainedNetPath = self.getTrainedNetPath('', true);
            trainedNetDir = fileparts(trainedNetPath);
            if ~exist(trainedNetDir, 'dir')
                mkdir(trainedNetDir);
            end
            save(trainedNetPath, 'net', 'tr', 'accuracy');
        end
        
        function loadNet(self)
            if ~self.testOnly
                trainedNetFile = self.getTrainedNetPath('', true);
            else
                trainedNetFile = self.existNetPath;
            end
            S = load(trainedNetFile);
            self.net = S.net;
            self.tr = S.tr;
            self.accuracy = S.accuracy;
        end
        
        function out = eval(self, imset, ids, overwrite)
            if nargin < 4
                overwrite = false;
            end
            
            resultsPath = self.getResultsPath(imset);
            if ~exist(resultsPath, 'dir')
                mkdir(resultsPath);
            end
            
            if isempty(self.net)
                self.loadNet();
            end
            
            opts.spix = true;
            
            % Here it crashes 
            [~, out] = self.dh.iterIds(imset, ids, opts, @evalSingleIm);

            function out = evalSingleIm(in, ~, getStruct)
                out.postp = 0;
                if getStruct
                    return;
                end
                
                if ~overwrite && ~self.OVERWRITE && ...
                        exist(self.getResultsPath(imset, in.id), 'file')
                    return;
                end

                featureTic = tic;

                features = self.loadFeatures(imset, in.id);
                featureVects = cell2mat(struct2cell(features));
                spixIds = fieldnames(features);
                spixNums = SuperPixUtils.getNumsFromIds(spixIds);
                spixSubset = spixNums;

                confs = self.net(featureVects');
                spixLabel = false(in.numSupers, 1);
                imSize = size(in.superPix);
                imLabel = false(imSize);
                imConf = zeros(imSize);
                
                spixUtils = SuperPixUtils(self.dh, in.imset, in.id);
                spixUtilsSubset = spixUtils.getSpixSubset();
                assert(isequal(spixSubset, spixUtilsSubset), ...
                    'Feature vectors do not match spix subset');
                
                spixUtils.iterSpix(in, [], @handleSpix);
                function handleSpix(~, ~, spin)
                    conf = confs(spin.i);
                    
                    label = conf > self.THRESH;
                    imLabel(spin.spixIndex) = label;
                    spixLabel(spin.i) = label;
                    
                    imConf(spin.spixIndex) = conf;
                end
                labeledSupers = spixSubset(spixLabel);
                
                featureTime = toc(featureTic);
                fprintf('Feature time for: %s = %f\n', in.id, featureTime);

                postpTic = tic;
                in.spixLabel = spixLabel;
                in.labeledSupers = labeledSupers;
                [imLabel, imConf, rest, postpRet] = ...
                    self.postprocess(imLabel, imConf, in, featureVects);
                out.postp = postpRet;
                postpTime = toc(postpTic);
                fprintf('Postp time for: %s = %f\n', in.id, postpTime);
                self.saveResults(imset, in.id, imLabel, imConf, rest);
            end            
        end
        
        function [imLabel, imConf, rest, ret] = ...
                postprocess(self, imLabel, imConf, in, featureVects) %#ok<*INUSD,INUSL>
            rest = struct;
            ret = struct;
        end
        
        function saveResults(self, imset, id, imLabel, imConf, rest)
            if nargin < 6
                rest = struct;
            end
            rest.imLabel = imLabel;
            rest.imConf = imConf; %#ok<STRNU>
            resultsPath = self.getResultsPath(imset, id);
            
            disp('generic netline 593');
            disp(resultsPath);
            
            save(resultsPath, '-struct', 'rest');
        end
        
        function res = loadResults(self, imset, id)
            resultsPath = self.getResultsPath(imset, id);
            res = load(resultsPath);
        end
        
        function prepAndEval(self, varargin)
            %self.prepFeatures(varargin{:});
            self.eval(varargin{:});
        end
        
        function [out, time] = iterResults(self, imset, ids, func, opts)
            if nargin < 5
                opts.noIm = true;
                opts.noPos = true;
            end
            
            opts.cellOnly = true;
            [out, ~, time] = self.dh.iterIds(imset, ids, opts, @handleSingleIm);
            function ret = handleSingleIm(in, ~, ~)   
                res = self.loadResults(in.imset, in.id);
                ret = func(in, res);
            end
        end
        
        function suffix = writeImConfs(self, imset, ids, labelsOnly, outDir, op)
            if nargin < 5
                outDir = '';
            end
            if nargin < 4
                labelsOnly = false;
            end

            suffix = '';
            outPath = self.getOutSubPath(imset, outDir, labelsOnly);
            if self.NESTED_OUT
                [outPath, suffix] = GenericNet.getUniquePath(outPath);
            end
            mkdir(outPath);

            numIds = numel(ids);
            for i = 1:numIds
                id = ids{i};

                res = self.loadResults(imset, id);
                if ischar(labelsOnly)
                    imLabel = res.(labelsOnly);
                elseif labelsOnly
                    imLabel = res.imLabel;
                else
                    imLabel = res.imConf;
                end
                
                if nargin >= 6
                    imLabel = op(imLabel);
                end
                
                outFile = fullfile(outPath, strcat(id, '.png'));
                imwrite(imLabel, outFile, 'png');
            end
        end
        
        function evalSeg(self, labelsOnly, suffix)
            if nargin < 3
                suffix = '';
            end
            imset = self.dh.SET_TRAIN;
            outPath = self.getOutSubPath(imset, '', labelsOnly);
            outPath = [outPath suffix];
            logPath = fullfile(outPath, 'eval.txt');
            utilitiesPath = fullfile(self.dh.PROJECT_PATH, 'code/devkit', ...
                'python', 'utilities.py'); 
            cmd = sprintf('python "%s" "%s"', utilitiesPath, outPath);
            system(cmd);
            logf = fopen(logPath, 'a');
            c = onCleanup(@() fclose(logf));
            exclude = {};
            changedProps = self.getChangedProps(true, exclude);
            changedPropStr = self.genPropStr(changedProps, {}, ...
                false, ' = ', '\n');
            allProps = self.getAllProps();
            propStr = self.genPropStr(allProps, {}, false, ' = ', '\n');
            dispVal = evalc('self');
            dhDisp = evalc('self.dh');
            trDisp = evalc('self.tr');
            fprintf(logf, '\n\nChanged:\n%s\n\nAll:\n%s\n%s\n%s\n%s', ...
                changedPropStr, propStr, dispVal, dhDisp, trDisp);
            suffix = '';
            if labelsOnly
                suffix = 'Bin';
            end
            hostname = evalc('!hostname');
            autogenSub = '';
            if ~isempty(strfind(hostname, 'cluster'))
                autogenSub = 'csl';
            end
            testProp = self.EVAL_TEST_PROP;
            if ~any(ismember(changedProps, testProp)) && isprop(self, testProp)
                logProps = cat(1, changedProps, testProp);
            else
                logProps = changedProps;
            end
            logCopyName = ['eval' suffix '-' self.genPropStr(logProps) '.txt'];
            logCopyDir = fullfile(self.dh.PROJECT_PATH, ...
                'code/bevscores/autogen', autogenSub, ...
                self.getNetSubdir(imset), testProp);
            mkdir(logCopyDir);
            logCopyPath = fullfile(logCopyDir, logCopyName);
            copyfile(logPath, logCopyPath);
        end
        
        function runAll(self, overwrite, prepFeat, train, eval, write, trainIds, testIds, testImset)
            allTic = tic;
            if nargin < 2
                overwrite = false;
            end
            if nargin < 7
                trainIds = self.trainIds;
            end
            if nargin < 8
                testIds = self.testIds;
            end
            if nargin < 9
                numTrainIms = numel(self.dh.getIds(self.dh.SET_TRAIN));
                % -5 since maybe excluding a couple on purpose
                if numel(trainIds) < numTrainIms - 5
                    testImset = self.dh.SET_TRAIN;
                else
                    testImset = self.dh.SET_TEST;
                end
            end
            if nargin < 3 || prepFeat
                self.prepFeatures(self.dh.SET_TRAIN, trainIds, overwrite);
                self.prepFeatures(testImset, testIds, overwrite);
                self.prepTrainData(overwrite);
            end
            if nargin < 4 || train
                self.train(overwrite);
            end
            if nargin < 5 || eval
                self.eval(testImset, testIds, overwrite);
            end
%             if nargin < 6 || write
            if nargin > 5 && write
                suffix = self.writeImConfs(testImset, testIds, true);
                suffix2 = self.writeImConfs(testImset, testIds, false);
                if strcmp(testImset, self.dh.SET_TRAIN)
                    self.evalSeg(true, suffix);
                    self.evalSeg(false, suffix2);
                end
            end
            allTime = toc(allTic);
            fprintf('Total running time: %f\n', allTime);
        end
    end
    
    methods(Static)
        function [uniquePath, suffix] = getUniquePath(path)
            basePath = path;
            suffix = '';
            i = 2;
            while exist(path, 'file')
                suffix = ['-' num2str(i)];
                path = [basePath suffix];
                i = i + 1;
            end
            uniquePath = path;
        end
    end
end
