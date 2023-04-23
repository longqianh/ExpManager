
function avi2gif(vidpath,delayTime)
    if nargin<2
        delayTime=0.05;
    end

    % Load the video file
    video = VideoReader(vidpath);
    
    % Extract each frame from the video and convert it to grayscale
    for i = 1:video.NumFrames
        frame = readFrame(video);
%         grayframe = rgb2gray(frame);   
        % Create a plot of the grayscale frame without axis
        imshow(frame); % plot the grayscale image
        colormap(gray); % set the colormap
        axis off; % remove the axis
        % Convert the plot to a GIF image
        if i == 1
            % Create the GIF file and write the first frame
            set(gcf, 'Color', 'w','MenuBar','none','ToolBar','none','resize','off'); % set the figure background to white
            set(gca, 'Position', [0 0 1 1]); % set the axis position to cover the entire figure
            set(gca,'visible','off')

            imwrite(rgb2gray(frame2im(getframe(gcf))), 'output.gif', 'DelayTime', delayTime, 'LoopCount', inf);
        else
            % Append subsequent frames to the existing GIF file
            imwrite(rgb2gray(frame2im(getframe(gcf))), 'output.gif', 'DelayTime', delayTime, 'WriteMode', 'append');
        end
        
        hold off; % release the plot
    end
end