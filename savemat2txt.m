
% 
% Copyright (C) 2015  Xiaozhi Chen, Kaustav Kundu, Yukun Zhu, Andrew Berneshawi, Huimin Ma, Sanja Fidler, Raquel Urtasun
% Website: http://www.cs.toronto.edu/objprop3d/
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


function savemat2txt(mat, fname, ftype, dtype)
%SAVEMAT2TXT Summary of this function goes here
%   Detailed explanation goes here

assert(ndims(mat) == 2, 'The matrix should contain 2 dimensions');

if nargin < 3
    ftype = 'txt';
end

[rows, cols] = size(mat);


fid = fopen(fname, 'w');
if strcmp(ftype, 'txt')
    fprintf(fid, '# Matrix\n');
    fprintf(fid, 'WIDTH %d\nHEIGHT %d\n', cols, rows);
    for i = 1 : rows
        fprintf(fid, '%d ', mat(i,:));
        fprintf(fid, '\n');
    end
elseif strcmp(ftype, 'bin')
    if nargin < 4
        dtype = 'float';
    end
    fwrite(fid, [rows; cols], dtype);
    fwrite(fid, mat, dtype);
else
    error('wrong file type\n');
end

fclose(fid);


end

