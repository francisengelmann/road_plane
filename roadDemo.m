
% 
% Copyright (C) 2015  Xiaozhi Chen, Kaustav Kundu, Yukun Zhu, Andrew Berneshawi, Huimin Ma, Sanja Fidler, Raquel Urtasun
% Website: http://www.cs.toronto.edu/objprop3d/
% 
% with modifications by Francis Engelmann
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

mfilepath=fileparts(which(mfilename));
addpath(fullfile(mfilepath, 'devkit'))
addpath(fullfile(mfilepath, 'matlab'))

% NOTE that you need spsstereo, found here:
% http://ttic.uchicago.edu/~dmcallester/SPS/
% Download and compile it, and put an "spsstereo" folder in the ../external
% folder, with the binary called "spsstereo" in the root of that folder
% You also need the neural network toolbox of Matlab

%% Testing only
% We provide a model trained on KITTI road training set, the trained net
% is in ./data/trainedNet.mat
% If you just want to evaluate, lets say on KITTI object training set,
% you just need to create a soft link of the images in
% data/object_training/, so that it has left and right images in image_2 
% and image_3 respectively. Then you can run:

dataDir = 'data'; outDir = 'data';
dh = DataHandler(outDir, dataDir);
roadNet = RoadNet(dh, false, false, 'data/trainedNet.mat');

% Remove frames for which plane computation fails 
% roadNet.testIds{1}=[] % Example: this would remove frame 000000
% testIds_new = roadNet.testIds(~cellfun(@isempty, roadNet.testIds));

roadNet.testIds = testIds_new;
roadNet.prepAndEval(dh.SET_TEST, roadNet.testIds);

% This will output the road planes of the evaluated dataset to:
% data/object_training/roadnet/road/eval/
% After this, you can read the road planes in that folder using function
% computeRoadPlane.m


%% Training and testing
% We use the KITTI road dataset to train a road classifier and test on
% object dataset. It will default to looking for data under ./data
% Format should be as if KITTI data was extracted into ./data folder, so
% data/road_training and data/object_training

%dh = DataHandler();
%roadNet = RoadNet(dh);
%roadNet.runAll()

% Note that the first run will be pretty slow as it has to generate all the
% stereo files, but it will be faster after that if you tweak something and
% rerun.
% The DataHandler class is actually very configurable and you can change
% all sorts of paths, almost nothing is hardcoded so take a look at the
% fields and set them as you like
