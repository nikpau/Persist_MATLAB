classdef Viz
    %VIZ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        f1
        f2
        f3
        f4
        plots
        plots2
        plots3
        plots4
        numShips
        shipToViz = 5;
        zoom = 2000;
%         Frames % for video record
    end
    
    methods
        function o = Viz(numShips,river,ships,SimSet)
            vizutm = SimSet.vizInUTM;
            o.numShips = numShips;
%             o.Frames = [];
            
            close all
            clf
            
            if vizutm
                o.f1 = subplot(4,1,1);
                o.f2 = subplot(4,1,2);
                o.f3 = subplot(4,1,3);
                o.f4 = subplot(4,1,4);

            else
                o.f1 = subplot(4,1,1);
                o.f2 = subplot(4,1,2);
                o.f3 = subplot(4,1,3);
                o.f4 = subplot(4,1,4);
            end
            x0=10;
            y0=100;
            width=1500;
            height=800;
            set(gcf,'position',[x0,y0,width,height])
            set(gcf, 'Position', get(0, 'Screensize'));
            
            if vizutm
                
                o.f1.Position = [0.02 0.08 0.5 0.91];
                o.f2.Position = [0.55 0.7 0.42 0.25];
                o.f3.Position= [0.55 0.4 0.42 0.25];
                o.f4.Position= [0.55 0.1 0.42 0.25];
            end

            %     ButtonH=uicontrol('Parent',figure,'Style','pushbutton','String','View Data','Units','normalized','Position',[0.0 0.5 0.4 0.2],'Visible','on');
            subplot(o.f1)
            % RIVER
            [X,Y] = meshgrid(1:20:20*length(river.PointCoordinates),0:20:500);
            if vizutm
%                 surf(river.PointCoordinates(:,:,2),river.PointCoordinates(:,:,1),river.waterDepth + 10)
                surf(river.PointCoordinates(:,:,2),river.PointCoordinates(:,:,1),river.streamVel + 10)
                view(90, -90)
            else
                surf(X,Y,river.waterDepth + 10)
                view(0, -90)
            end
            shading interp
            hold on
           
           if vizutm
           pbaspect([ 1 1 1])
           end
            % all other ships
            for i=1:numShips
                
%                 plots(i) = scatter(0,0,50,'s','filled','DisplayName',num2str(i));
%                plots(i) =  rectangle('Position',[0 0 0 0]);
                if vizutm                  
                    plots(i) = plot(o.transformBoxInUTM(ships.heading.box{i},river));
                else
                    plots(i) = plot(ships.heading.box{i});
                end
              
                hold on
            end
%             legend show
            
            % river
            if vizutm
                xy = river.getPosUTMfromPos([10000;250]);
                plots(numShips + 1) = plot(xy(1),xy(2),'*');
                plots(numShips + 2) = plot(xy(1),xy(2),'y*');
                
%                 xlim([xy(1) - 5000 xy(1) + 5000])
%                 xlim([xy(2) - 5000 xy(2) + 5000])
               
            else
                
                plots(numShips + 1) = plot(0,0,'*-');
                plots(numShips + 2) = plot(0,0,'*-');
                
                xlim(ships.x(o.shipToViz)+[-2000 2000])
                ylim([0 520])
            end
            
            subplot(o.f2)
            
            % VIZ for LatConPol Observations
            plots2(1) = scatter(0,0,50,'s','filled','DisplayName',num2str(1));
            hold on
            plots2(2) = scatter(0,0,50,'s','filled','DisplayName',num2str(2));
            plots2(3) = scatter(0,0,50,'s','filled','DisplayName',num2str(2));
            xlim([-200 1000])
            ylim([-500 500])

             % VIZ for LonConPol Observations   
%             for i=1:numShips                
%                plots2(i) =  plot(0,0,'DisplayName','Engine Power');              
%                hold on
%             end
%             legend show
            
            
            
            subplot(o.f3)
            for i=1:numShips                
               plots3(i) =  plot(0,0,'DisplayName','Vessel velocity in m/s');              
               hold on
            end
            title('Vessel velocity in m/s')
            
           subplot(o.f4)
            for i=1:numShips                
               plots4(i) =  plot(0,0,'DisplayName','Water depth below keel in m');              
               hold on
            end
            title('Water depth below keel in m')

%             legend show
            
            o.plots = plots;
            o.plots2 = plots2;
            o.plots3 = plots3;
            o.plots4 = plots4;
            
            

        end
        
        
        %% VISUALIZE 
        function [o, SimSet] = visualize(o,ships,river,SimSet)

          vizutm = SimSet.vizInUTM;
          
            for i=1:o.numShips
                if vizutm
                    o.plots(i).Shape = o.transformBoxInUTM(ships.heading.box{i},river);
                else
                    o.plots(i).Shape = ships.heading.box{i};
                    %             o.f1.XLim = [ships.x(1)-4000  ships.x(1)+ 500];
                end
            end
            
            if vizutm
               
                   o.f1.XLim = [ships.x_utm(o.shipToViz)-o.zoom  ships.x_utm(o.shipToViz)+ o.zoom];
                   o.f1.YLim = [ships.y_utm(o.shipToViz)-o.zoom  ships.y_utm(o.shipToViz)+ o.zoom];
            else
                  o.f1.XLim = ships.x(o.shipToViz)+[-1800 1800];
                o.f1.YLim = [0 520];
            end
                        vizutm = SimSet.vizInUTM;

          
            
            
            
            
            
            
            % river
            if vizutm
                try
                for i=1:length(ships.latConPol(o.shipToViz).basePointsL(:,1))  
                    basePointsL_utm(:,i) = river.getPosUTMfromPos(ships.latConPol(o.shipToViz).basePointsL(i,:));
                    basePointsU_utm(:,i) = river.getPosUTMfromPos(ships.latConPol(o.shipToViz).basePointsU(i,:));          
                end
                
                
                o.plots(o.numShips + 1).XData = basePointsL_utm(1,:);
                o.plots(o.numShips + 1).YData = basePointsL_utm(2,:);
                o.plots(o.numShips + 2).XData =basePointsU_utm(1,:);
                o.plots(o.numShips + 2).YData = basePointsU_utm(2,:);
                end
            else
                o.plots(o.numShips + 1).XData = ships.latConPol(o.shipToViz).basePointsL(:,1);
                o.plots(o.numShips + 1).YData = ships.latConPol(o.shipToViz).basePointsL(:,2);
                o.plots(o.numShips + 2).XData = ships.latConPol(o.shipToViz).basePointsU(:,1);
                o.plots(o.numShips + 2).YData = ships.latConPol(o.shipToViz).basePointsU(:,2);
            end
            
           % VIZ for LatConPol Observations           
            o.plots2(1).XData = [ships.latConPol(o.shipToViz).TTCvectorUp];
            o.plots2(1).YData = [ships.latConPol(o.shipToViz).POSvectorUp];
            o.plots2(2).XData = [ships.latConPol(o.shipToViz).TTCvectorLow];
            o.plots2(2).YData = [ships.latConPol(o.shipToViz).POSvectorLow];

           % VIZ for LonConPol Observations          
            
%             for i=1:o.numShips
%                 o.plots2(i).XData = [o.plots2(i).XData SimSet.t]; 
%                 o.plots2(i).YData = [o.plots2(i).YData ships.ax(i)]; 
%                 
%             end

            for i=1:o.numShips
                o.plots3(i).XData = [o.plots3(i).XData SimSet.t]; 
%                 o.plots3(i).YData = [o.plots3(i).YData ships.lonConPol(i).distShipAhead]; 
                 o.plots3(i).YData = [o.plots3(i).YData abs(ships.vx(i))]; 
%                 o.plots3(i).YData = [o.plots3(i).YData ships.ax(i)]; 
%                 o.plots3(i).YData = [o.plots3(i).YData ships.P(i)]; 
%                 o.plots3(i).YData = [o.plots3(i).YData ships.ssq(i)]; 
%                 o.plots3(i).YData = [o.plots3(i).YData ships.heading.angle(i)]; 
%                 o.plots3(i).YData = [o.plots3(i).YData ships.y(i)]; 
                
            end
            
           for i=1:o.numShips
                o.plots4(i).XData = [o.plots4(i).XData SimSet.t]; 
%                 o.plots4(i).YData = [o.plots4(i).YData ships.lonConPol(i).distShipAhead]; 
%                 o.plots4(i).YData = [o.plots4(i).YData ships.vx(i)]; 
%                 o.plots4(i).YData = [o.plots4(i).YData ships.ax(i)]; 
%                  o.plots4(i).YData = [o.plots4(i).YData ships.P(i)]; 
                % water depth under keel
                o.plots4(i).YData = [o.plots4(i).YData ships.ssq(i)-river.waterDepth(round(ships.y(i)/20),round(ships.x(i)/20))]; 
%                 o.plots4(i).YData = [o.plots4(i).YData ships.heading.angle(i)]; 
%                 o.plots4(i).YData = [o.plots4(i).YData ships.y(i)]; 
                
           end
%         o.Frames = [o.Frames getframe(gcf)];
            
         pause(0.002)
         if ~ishghandle(o.f1)
             SimSet.isDone = 1;
         end
        end
        
        function box = transformBoxInUTM(~,box,river)
            newVertices = zeros(4,2);
            for i=1:4
                newVertices(i,:) = river.getPosUTMfromPos(box.Vertices(i,:))';
            end
            box = polyshape(newVertices(:,1)', newVertices(:,2)');
                
        end
    end
end

