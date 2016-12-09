### Road plane computation

Reassambled bundle for road-plane computation from [1] for my conveniance.

### Usage
- Create a symbolic link to your dataset in ```./data``` named data_object (see ```./matlab/DataHandler.m```).
- In ```./matlab/DataHandler.m``` appl the following pathes:
    - line 13: do you have subdirs?
    - line 15: is there a calibration file per sequence?
    - line 18: this is the name of you symbolic link in ```./data```.
    - line 239: is this path correct?
- You need sps-stereo (disparity estimation method) to be installed under ```../spsstereo```.
- Run convertMAT2TXT.m to extract the plane parameters into a human-readable txt.

### References

[1] Xiaozhi Chen, Kaustav Kundu, Yukun Zhu, Andrew Berneshawi, Huimin Ma, Sanja Fidler, Raquel Urtasun. 
3D Object Proposals for Accurate Object Class Detection 
Neural Information Processing Systems (NIPS), Montreal, Canada, 2015 
