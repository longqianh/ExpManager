classdef MyClass < handle
    properties
        canvas
        canvas_t
    end
    
    methods
        function obj = MyClass()
            obj.canvas = figure('Color','White','Name',"Propagation Figures","Visible","off");
            obj.canvas_t = tiledlayout(obj.canvas,'flow','TileSpacing','none','Padding','none');
        end
        
        function plotData(obj, data)
            % Plot the data on the tiled layout
            nexttile(obj.canvas_t)
            plot(data)
            obj.canvas.Visible=1;
        end
    end
end
