classdef Image_Analysis_Horizontal_Hong < handle
    % Written by Moo Sun Hong on Nov 6, 2018
    % Editted by Moo Sun Hong on May 27, 2021
    % Hong et al., A Droplet-based Evaporative System for the Estimation of Protein Crystallization Kinetics
    % Detect boundary of droplet from horizontal view images
    % Calculate volume based on symmtery assumption
    properties
        filt=50;
        threshold=-0.05;
        xclear=10;
        yclear=3;
        yfocus=0;
    end
    methods
        function volume = processframe(obj,frame)
            figure(1); imshow(frame); % Figure 2a
            
            lowfreq = imgaussfilt(frame,obj.filt);
            highfreq = frame-lowfreq;
            framethred = highfreq<obj.threshold;
            figure(2); imshow(framethred); % Figure 2b
            
            framethred = imopen(framethred,strel('rectangle',[obj.yclear;obj.xclear]));
            framethred = imclose(framethred,strel('rectangle',[obj.yclear;obj.xclear]));
            framethred = imfill(framethred,'holes');
            figure(3); imshow(framethred); % Figure 2c
            
            stats = regionprops(framethred,'BoundingBox','ConvexArea','ConvexImage');
            [~,ind] = max([stats.ConvexArea]);
            figure(4); imshow(frame); % Figure 2d
            if ~isempty(ind)
                area = stats(ind).ConvexImage;
                rbound = zeros(size(area,1),2);
                for ii = 1:size(area,1)
                    rbound(ii,1) = find(area(ii,:),1);
                    rbound(ii,2) = find(area(ii,:),1,'last');
                end
                convind1 = find(rbound(:,1) == min(rbound(:,1)),1);
                convind2 = find(rbound(:,2) == max(rbound(:,2)),1);
                convind = max(convind1,convind2);
                yup = stats(ind).BoundingBox(2)-0.5+convind;
                convind = yup-stats(ind).BoundingBox(2)+0.5;
                hold on;
%                 plot(stats(ind).BoundingBox(1)-1+[rbound(1:convind,1);flip(rbound(1:convind,2));rbound(1,1)]...
%                     ,stats(ind).BoundingBox(2)-0.5+[(1:convind),flip((1:convind)),1]','b'...
%                     ,(1:size(frame,2))',obj.yfocus*ones(size(frame,2),1),'g','LineWidth',2);
                rbound = rbound(convind:end,:);
                plot(stats(ind).BoundingBox(1)-1+[rbound(:,1);flip(rbound(:,2));rbound(1,1)]...
                    ,stats(ind).BoundingBox(2)-1.5+convind+[(1:length(rbound)),flip((1:length(rbound))),1]','r','LineWidth',2);
                hold off;
                volume = pi/4*sum((rbound(:,2)-rbound(:,1)+1).^2); ... % Volume [px^2]
            else
                volume = nan;
            end            
        end
    end
end