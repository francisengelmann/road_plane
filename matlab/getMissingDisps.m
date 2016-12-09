function missDisps = getMissingDisps(ids, imset)

if nargin < 2
    imset = 'test';
end

numIds = size(ids, 1);
sids = strcat('id', ids);
dispsCell = cell(numIds, 1);
parfor i = 1:numIds
    id = ids{i};
    disp = getDisparity(imset, id);
    missDisp = disp == 0;
    dispsCell{i} = missDisp;
end

missDisps = cell2struct(dispsCell, sids);

end
