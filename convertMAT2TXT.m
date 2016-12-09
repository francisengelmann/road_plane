dir_in = '/media/engelmann/6de91958-d0ea-4752-90ec-47c2b0046cce/work/francis/kitti/data_tracking/training/0019/roadnet/road/eval';
dir_out = '/media/engelmann/6de91958-d0ea-4752-90ec-47c2b0046cce/work/francis/kitti/data_tracking/training/0019/planes';

first_frame = 0
last_frame = 5000

for id  = first_frame:last_frame
    try 
        load(sprintf('%s/%06d.mat', dir_in, id));
        plane_file=sprintf('%s/%06d.txt', dir_out, id);
        if plane.normal(2) > 0
            % make sure the normal is always facing up
            plane.normal = -plane.normal;
        end
        plane = [plane.normal, -plane.point * plane.normal'];
        plane = plane ./ norm(plane(1:3));
        savemat2txt(plane(:)', plane_file);
    catch
        error('Stopped at frame %d',id);
    end
end