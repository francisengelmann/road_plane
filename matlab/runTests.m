function runTests(trainIds, testIds, testNo, argVals)

if nargin < 4
    argVals = [];
end
in.trainIds = trainIds;
in.testIds = testIds;
in.argVals = argVals;

%imset = 'train';

if any(testNo == 1)
    constVals = [30 50 100 200 300 500 700 1000 1e4 1e5 1e6];
    prop = 'dh.depthLims.lim';
    runTest(in, constVals, prop);
end

if any(testNo == 2)
    constVals = {{1000 1000 1000} {50 20 1000} {50 20 1e6} {100 50 1e6}};
    prop = 'dh.depthLims.postLims';
    runTest(in, constVals, prop)
end

if any(testNo == 3)
    constVals = [10 20 50 100 200 300 400 500 750 1000];
    prop = 'maxFail';
    runTest(in, constVals, prop);
end

if any(testNo == 4)
    constVals = [250 500 750 1000 1500 2000 3000 3500];
    prop = 'dh.numSpix';
    runTest(in, constVals, prop);
end

if any(testNo == 5)
    constVals = [500 1000 1500 2000 3000];
    vals = getVals(argVals, constVals);
    posPercents = [0.3 0.6 0.7 0.8 0.9 0.99];
    for i = 1:numel(vals)
        for j = 1:numel(posPercents)
            dh = DataHandler();
            dh.numSpix = vals(i);
            roadNet = RoadNet(dh, trainIds, testIds);
            roadNet.MIN_POS_PERCENT = posPercents(j);
            roadNet.EVAL_TEST_PROP = 'MIN_POS_PERCENT';
            roadNet.runAll();
            killPool();
        end
    end
end

if any(testNo == 6)
    constVals = {true, false, true};
    props = {'fitPlane', 'fitPlane', 'bwConn'};
    propFolder = 'postp';
    runTest(in, constVals, props, @(varargin)false, propFolder);
end

if any(testNo == 7)
    constVals = [1 2 3 5 10 20 50];
    layers = {[], [30], [30 30], [50], [50 100], [20 30 20]};
    vals = getVals(argVals, constVals);
    for i = 1:numel(vals)
        for j = 1:numel(layers)
            dh = DataHandler();
            roadNet = RoadNet(dh, trainIds, testIds);
            roadNet.EVAL_TEST_PROP = 'contextSize';
            roadNet.useContextFeat = true;
            roadNet.contextSize = vals(i);
            roadNet.setTotalFeatures();
            layer = layers{j};
            if ~isempty(layer)
                roadNet.layers = layer;
            end
            roadNet.runAll();
            killPool();
        end
    end
end

end

function runAll(roadNet, varargin)
    roadNet.runAll();
end

function runTest(in, constVals, prop, func, propFolder, varargin)
    if nargin < 4
        func = @runAll;
    end
    if nargin < 5
        propFolder = prop;
    end
    
    function roadNet = setup(~, ~, ~, in)
        dh = DataHandler();
        roadNet = RoadNet(dh, in.trainIds, in.testIds);
        roadNet.EVAL_TEST_PROP = propFolder;
    end
    
    tester = Tester(@setup, Tester.nop, true);
    tester.runTest(in, constVals, prop, func, varargin{:})
end

function vals = getVals(varargin)
    vals = Tester.getVals(varargin{:});
end

function killPool()
    Tester.killPool();
end
