classdef Image_Analysis_Vertical_Hong < handle
    % Written by Moo Sun Hong on Jul 13, 2020
    % Editted by Moo Sun Hong on May 27, 2021
    % Hong et al., A Droplet-based Evaporative System for the Estimation of Protein Crystallization Kinetics
    % Detect boundary of crystals from vertical view images
    properties
        filt=250;
        threshold=0.01;
        clear=3;
    end
    methods
        function area = processframe(obj,frame)
           figure(1); imshow(frame); % Figure 3a
           
           lowfreq = imgaussfilt(frame,obj.filt);
           highfreq = frame-lowfreq;
           framethred = highfreq<obj.threshold;
           figure(2); imshow(framethred); % Figure 3b
             
            stats = regionprops(framethred,'PixelList'); pixel = {stats.PixelList};
            for ii = 1:length(pixel)
                if sum(sqrt(sum((pixel{ii}-800).^2,2))>660) > 0
                    framethred((pixel{ii}(:,1)-1)*1600+pixel{ii}(:,2)) = 0;
                end
            end
            figure(3); imshow(framethred); % Figure 3c
            
            framethred = imdilate(framethred,strel('disk',obj.clear));
            framethred = imfill(framethred,'holes');
            framethred = imerode(framethred,strel('disk',obj.clear));
            figure(4); imshow(framethred); % Figure 3d
            
            framethred = imopen(framethred,strel('disk',obj.clear));
            figure(5); imshow(framethred); % Figure 3e
            
            stats = regionprops(framethred,'Area');
            area = [stats.Area];            
            boundary = bwboundaries(framethred,'noholes');
            figure(6); imshow(frame); % Figure 3f
            hold on; cellfun(@(x)plot(x(:,2),x(:,1),'r','LineWidth',2),boundary); hold off;
        end
    end
end