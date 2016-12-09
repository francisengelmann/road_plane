classdef Tester < matlab.mixin.Copyable

    properties(Constant)
        nop = @(x, varargin) x;
    end
    
    properties
        setupPtr;
        tearDownPtr;
        killPoolFlag;
    end
    
    methods
        function tester = Tester(setupFunc, tearDownFunc, killPoolFlag)
            self = tester;
            nopf = self.nop;
            if nargin < 1 || isempty(setupFunc)
                setupFunc = nopf;
            end
            if nargin < 2
                tearDownFunc = nopf;
            end
            if nargin < 3
                killPoolFlag = false;
            end
            self.setupPtr = setupFunc;
            self.tearDownPtr = tearDownFunc;
            self.killPoolFlag = killPoolFlag;
        end
        
        function obj = setup(self, obj, prop, val, in)
            obj = self.setupPtr(obj, prop, val, in);
        end
        
        function tearDown(self, res)
            self.tearDownPtr(res);
        end
        
        function res = runTest(self, in, constVals, props, func, obj)
            vals = Tester.getVals(in.argVals, constVals);
            if ~iscell(vals)
                vals = num2cell(vals);
            end
            if ~iscell(props)
                props = repmat({props}, 1, numel(vals));
            end
            if nargin < 6
                obj = [];
            end
            for i = 1:numel(vals)
                prop = props{i};
                val = vals{i};
                obj = self.setup(obj, prop, val, in);
                setCmd = ['obj.' prop ' = val;'];
                eval(setCmd);
                if nargout(func)
                    res = func(obj, prop, val, in);
                else
                    func(obj, prop, val, in);
                    res = false;
                end
                self.tearDown(res);
                if self.killPoolFlag
                    Tester.killPool();
                end
            end
        end
    end
    
    methods(Static)
        function vals = getVals(argVals, constVals)
            if isempty(argVals)
                vals = constVals;
            else
                vals = argVals;
            end
        end
        
        function killPool()
            delete(gcp('nocreate'));
            if isunix()
                system('pkill -u bernesha -f "MATLAB -dmlworker"');
            end
        end
    end
end

