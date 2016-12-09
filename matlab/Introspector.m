classdef (Abstract) Introspector < matlab.mixin.Copyable

    properties
        defaultObj;
        disableInspect = false;
    end
    
    methods (Abstract)
        obj = constructDefault(self);
    end
    
    methods (Access = private)
        function excludeMap = getExcludeMap(self, exclude)
            if nargin < 2
                exclude = {};
            end
            exclude = cat(2, exclude, {'defaultObj'}, self.getExcludeProps());
            excludeMap = containers.Map(exclude, true(size(exclude)));
        end
    end
    
    methods
        function defaultObj = getDefaultObj(self)
            if self.disableInspect
                defaultObj = self;
            else
                if isempty(self.defaultObj)
                    self.defaultObj = self.constructDefault();
                end
                defaultObj = self.defaultObj;
            end
        end
        
        function exclude = getExcludeProps(~)
            % Exclude from both getChangedProps and getPropStr
            exclude = {};
        end
        
        function exclude = getExcludeChangedProps(~)
            % Exclude from only getChangedProps
            exclude = {};
        end
        
        function include = getIncludeProps(~)
            include = {};
        end
        
        function props = getAllProps(self)
            props = cat(1, fieldnames(self), self.getIncludeProps()');
        end
        
        function [changed, inds] = getChangedProps(self, props, exclude)
            if nargin < 2 || isequal(props, true)
                props = self.getAllProps();
            end
            if nargin < 3
                exclude = {};
            end
            exclude = cat(2, exclude, self.getExcludeChangedProps());
            excludeMap = self.getExcludeMap(exclude);
            defltObj = self.getDefaultObj();
            changed = {};
            inds = [];
            for i = 1:numel(props)
                prop = props{i};
                startSize = numel(changed);
                if ~isKey(excludeMap, prop)
                    val = self.getPropVal(prop);
                    default = defltObj.getPropVal(prop);
                    if isa(val, 'Introspector')
                        subChanged = val.getChangedProps();
                        subChanged = strcat([prop '.'], subChanged);
                        changed = cat(1, changed, subChanged);
                    elseif ~isequal(val, default)
                        if isstruct(val)
                            subChanged = self.getStructChanges(prop, ...
                                val, default);
                            changed = cat(1, changed, subChanged);
                        else
                            changed = cat(1, changed, prop);
                        end
                    end
                    if numel(changed) > startSize
                        inds = cat(1, inds, i);
                    end
                end
            end
            changed = unique(changed);
        end
        
        function changed = getStructChanges(~, prop, valStruct, defaultStruct)
            fields = fieldnames(valStruct);
            changed = {};
            for i = 1:numel(fields)
                field = fields{i};
                val = valStruct.(field);
                defHasField = isfield(defaultStruct, field);
                if defHasField
                    default = defaultStruct.(field);
                end
                if ~defHasField || ~isequal(val, default)
                    changed = cat(1, changed, field);
                end
            end
            changed = strcat([prop '.'], changed);
        end
        
        function val = getPropVal(self, prop)
            subProps = strsplit(prop, '.');
            val = self;
            for j = 1:numel(subProps)
                val = val.(subProps{j});
            end
        end
        
        function propStr = genPropStr(self, props, ops, shorten, assign, sep)
            if nargin < 3 || isempty(ops)
                ops = cell(numel(props), 1);
            end
            if nargin < 4
                shorten = true;
            end
            if nargin < 5
                assign = '=';
            end
            if nargin < 6
                sep = '+';
            end
            % Handle special chars
            sep = sprintf(sep);
            excludeMap = self.getExcludeMap();
            propStr = '';
            for i = 1:numel(props)
                prop = props{i};
                if ~isKey(excludeMap, prop)
                    val = self.getPropVal(prop);
                    op = ops{i};
                    if ~isempty(op)
                        val = op(val);
                    end
                    val = Introspector.val2str(val, shorten);
                    if shorten
                        propChars = Introspector.wordChars(prop);
                    else
                        propChars = prop;
                    end
                    propStr = [propStr, propChars, assign, val, sep]; %#ok<AGROW>
                end
            end
            propStr = propStr(1:end-1);
        end
        
        function propStr = genPropStrLines(self, varargin)
            % Convenience method
            defaults = {self.getAllProps(), {}, false, '=', '\n'};
            offset = numel(varargin);
            inputArgs = cat(2, varargin, defaults(offset+1:end));
            propStr = self.genPropStr(inputArgs{:});
        end
        
        function path = genPath(self, parent, prefix, props, fileId, ops, minProps)
            if nargin < 5
                fileId = '';
            end
            if nargin < 6
                ops = {};
            end
            if nargin < 7
                minProps = 1;
            end
            ext = '';
            if ~isempty(fileId)
                ext = '.mat';
            end
            if ~isempty(props)
                [changedProps, inds] = ...
                    self.getChangedProps(props(minProps+1:end));
                if ~isempty(ops)
                    % Select the subset of operations based on the changed
                    % props
                    inds = inds + minProps;
                    ops = cat(1, ops(1:minProps), ops(inds));
                end
                % Always include the first minProps but only include the rest if
                % they are changed
                props = cat(1, props(1:minProps), changedProps);
            end
            propStr = self.genPropStr(props, ops);
            if ~isempty(propStr) && ~isempty(prefix)
                % add hyphen
                prefix = [prefix '-'];
            end
            folder = [prefix propStr];
            path = fullfile(parent, folder);
            path = fullfile(path, [fileId ext]);
        end
    end
    
    methods(Static)
        function chars = wordChars(varname)
            varname = regexprep(varname, '([a-z])([A-Z])|(\w\.)(\w)', '$1_$2');
            toks = regexp(varname, '(?:^|_)(.)', 'tokens', 'ignorecase');
            % regexp for some reason nests each character in a cell
            toks = cellfun(@(x) cat(2, x{:}), toks, 'UniformOutput', false);
            if iscell(toks{1})
                chars = cellfun(@(x) strjoin(x, ''), toks, ...
                    'UniformOutput', false);
            else
                chars = strjoin(toks, '');
            end
            chars = lower(chars);
        end
        
        function valstr = val2str(val, shorten)
            if nargin < 2
                shorten = true;
            end
            try
                if isstruct(val)
                    fields = fieldnames(val);
                    val = struct2cell(val);
                    if ~shorten
                        vals = strtrim(cellstr(num2str(cell2mat(val))));
                        val = cat(2, fields, vals);
                        delims = repmat({':' '; '}, 1, numel(fields));
                        delims(end) = [];
                        val = strjoin(reshape(val', 1, numel(val)), delims);
                    end
                end
                if iscell(val)
                    if ischar(val{1})
                        if size(val, 1) > 1
                            val = val';
                        end
                        val = strjoin(val, ';');
                    else
                        val = cell2mat(val);
                    end
                end
                
                if ~ischar(val)
                    valstr = mat2str(val);
                    % Square brackets cause problems with python glob
                    % Change brackets to round brackets
                    valstr = strrep(strrep(valstr, '[', '('), ']', ')');
                else
                    valstr = val;
                end
            catch ME
                valstr = '';
            end
        end
    end
end

