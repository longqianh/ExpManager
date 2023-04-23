n=60;
imgs_path=cell(n,1);
for i=1:n
    imgs_path{i}=fullfile(ma.exp_save_dir,[num2str(i),'.bmp']);
end
%%

% create the videofile with 30 fps
% videofile = VideoWriter('output.avi','Grayscale AVI');
videofile = VideoWriter('output.mp4','MPEG-4');

videofile.FrameRate = 10;
% open the videofile
open(videofile);
% write the frames to the video
for i=1:length(imgs_path)
    % convert the image to a frame
    frame = imread(imgs_path{i});
    writeVideo(videofile, frame);
end
% close the videofile
close(videofile);