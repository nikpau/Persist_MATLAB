classdef River
    
    properties
        length
        basePointsL
        basePointsU
        minWaterUnderKeel
        
        streamVel
        waterDepth
        PointDataset
        PointCoordinates
    end
    
    methods
        function o = River()
            o.minWaterUnderKeel = 0.2;
            
            %load river data from inp file
            
            % river part from km 655 - 665
%             opts=detectImportOptions('MW_grid_500m_655_665_UTM32.inp','FileType','text');
%             data2 = readmatrix('MW_grid_500m_655_665_UTM32.inp',opts);

        % whole river
            opts=detectImportOptions('MW_grid_500m_655_852_UTM32.inp','FileType','text');
            data2 = readmatrix('MW_grid_500m_655_852_UTM32.inp',opts);
            
            index_ones = find(data2(:,1)==1);
            o.PointCoordinates = data2(index_ones(1):index_ones(2)-9,2:3); %Matrix with all point coordinates
            o.PointDataset = data2(index_ones(2):end,2:7); % Matrix with all point dataset (like stream velocity,...)
            
            o.streamVel = zeros(length(o.PointCoordinates),1);
            o.waterDepth = zeros(length(o.PointCoordinates),1);

            for i=1:length(o.PointCoordinates)
                o.streamVel(i) = sqrt(o.PointDataset(i,1)^2 + o.PointDataset(i,2)^2);
                o.waterDepth(i) = o.PointDataset(i,3);
            end
            
            o.streamVel  = reshape(o.streamVel,26,[]);
            o.waterDepth  = reshape(o.waterDepth,26,[]);            
            o.PointDataset = reshape(o.PointDataset,26,[],6);
            o.PointCoordinates = reshape(o.PointCoordinates,26,[],2);
            
            
            
            o.length = 20*length(o.PointCoordinates);
            
            %%%%%%%% CREATE BASEPOINTs %%%%%%%%
            
            % X VALUE rand distr with specific mean dist
            meanBasePointDist = 500;
            %             o.basePointsL = [0 sort(randperm(o.length,round(o.length/meanBasePointDist))) o.length];
            %             o.basePointsU = [0 sort(randperm(o.length,round(o.length/meanBasePointDist))) o.length];
            
            o.basePointsL = (0:meanBasePointDist:o.length);
            o.basePointsL = o.basePointsL + normrnd(0,10,[1 length(o.basePointsL)]);
            o.basePointsU = (0:meanBasePointDist:o.length);
            o.basePointsU = o.basePointsU + normrnd(0,10,[1 length(o.basePointsU)]);
            o.basePointsL(2,:) = 0;
            o.basePointsU(2,:) = 0;
            
            
            % Y VALUE with AR1 process
            AR1 = 20;
            for i=1:length(o.basePointsU)
                AR1 = 0.0100 + AR1 *  0.9995 + normrnd(0,sqrt( 0.3998));
                
                %                 if isber(i,o.basePointsL(1,:))
                %                     o.membasePointsL(2,find(o.basePointsL(1,:) == i)) = AR1 - 20 - randi([20 40]);
                temp = 55;
                o.basePointsL(2,i) = 250-randi(temp*[1 1]);
                %                 elseif ismember(i,o.basePointsU(1,:))
                %                     o.basePointsU(2,find(o.basePointsU(1,:) == i)) = AR1 - 20 + randi([20 40]);
                o.basePointsU(2,i) = 250+randi(temp*[1 1]);
                %                 end
            end
            
        end
        
        function   [xy_utm] = getPosUTMfromPos(o,xy)
            x = xy(1);
            y = xy(2);
            x_log = floor(x/20);
            xAlpha = x/20 - x_log;
            y_log = floor(y/20)+1;
            yAlpha = y/20 + 1 - y_log; 
            
            base_utm = o.PointCoordinates(y_log,x_log,:);
            
            % in x direction           
            diff = o.PointCoordinates(y_log,x_log + 1,:) - o.PointCoordinates(y_log,x_log,:);
            
            base_utm(:,:,1) = base_utm(:,:,1) + xAlpha * diff(:,:,1) ;
            base_utm(:,:,2) = base_utm(:,:,2) + xAlpha * diff(:,:,2) ;
            
            % in y direction           
            diff = o.PointCoordinates(y_log + 1,x_log,:) - o.PointCoordinates(y_log,x_log,:);
            x_utm = base_utm(:,:,1) + yAlpha * diff(:,:,1) ;
            y_utm = base_utm(:,:,2) + yAlpha * diff(:,:,2) ;
           
           xy_utm = [y_utm;x_utm];
            
                       
            
        end
        
        function  basePointsL = getClosestBasePointsL(o,x,dir,num)
            
            if dir == 1
                index = find(o.basePointsL(1,:) > x,1)-1;
            else
                index = find(o.basePointsL(1,:) > x,1);
            end
            isDone = 0;
            basePointsL = [];
            i = index;
            while ~isDone
                basePointsL = [basePointsL o.basePointsL(:,i)];
                i = i + dir;
                if abs(i - (index)) == num
                    isDone = 1;
                end
            end
        end
        
        function  basePointsU = getClosestBasePointsU(o,x,dir,num)
            
            if dir == 1
                index = find(o.basePointsU(1,:) > x,1)-1;
            else
                index = find(o.basePointsU(1,:) > x,1);
            end
            isDone = 0;
            basePointsU = [];
            i = index;
            while ~isDone
                basePointsU = [basePointsU o.basePointsU(:,i)];
                i = i + dir;
                if abs(i - (index)) == num
                    isDone = 1;
                end
            end
            
        end
        
        function  vStream = getMeanVStreamForVessel(o,ships, ID)
            
            minY_log = round((ships.y(ID)-ships.width_eff(ID)/2)/20);
            maxY_log = round((ships.y(ID)+ships.width_eff(ID)/2)/20);
            minX_log = round((ships.x(ID)-ships.length(ID)/2)/20)+10;
            maxX_log = round((ships.x(ID)+ships.length(ID)/2)/20)+10;
                        
            vStream = ships.dir(ID) *  mean(o.streamVel(minY_log:maxY_log,minX_log:maxX_log),'all');
%              vStream = 2;
        end
        
        function A = getRiverProfileForVessel(o,ships, ID)
            minX_log = round((ships.x(ID)-ships.length(ID)/2)/20)-10;
            maxX_log = round((ships.x(ID)+ships.length(ID)/2)/20)+10;
            
            n = maxX_log - minX_log;
            A = mean(trapz(o.waterDepth(:,minX_log:maxX_log))) * 20;

%             A = 1500;
        end

        function h = getWaterDepthForVessel(o,ships,ID)
            minX_log = round((ships.x(ID)-ships.length(ID)/2)/20)-10;
            maxX_log = round((ships.x(ID)+ships.length(ID)/2)/20)+10;
            minY_log = floor(ships.y(ID)/20);
            maxY_log = ceil(ships.y(ID)/20);
            
            h = mean(mean(o.waterDepth(minY_log:maxY_log,minX_log:maxX_log)));
%             h = 5;
        end
        function [lowerPoints, upperPoints] = detectFreeSpace(o,ID,ships,lookAheadDistance,numPoints)
            
            lookAheadDistance = lookAheadDistance/20;
            posX = round((ships.x(ID)-ships.length(ID) / 2)/20);
            ssq = ships.ssq(ID);
            
            lengthBlock = lookAheadDistance/numPoints;
            if ships.dir(ID) == 1
                start = round((floor(posX/lengthBlock)-1)*lengthBlock);
            else
                start = round((ceil(posX/lengthBlock)+1)*lengthBlock) - lookAheadDistance+1;
            end
            waterDepth = o.waterDepth(:,start:start + lookAheadDistance-1);              
            
            
            
            freeSpace = waterDepth -  (ssq + o.minWaterUnderKeel) < 0; % binary array ( 1 -> coast ; 0 -> freeSpace)
            
            % consider vessel length 
            for i=1:round(ships.length(ID)/20)
            freeSpace = [freeSpace zeros(26,1)] | [zeros(26,1) freeSpace];
            end
            
            lowerPoly = [];
            upperPoly = [];
            for i=1:size(freeSpace,2)
                [~,b]= max(diff(find(diff([1; freeSpace(:,i); 1])~=0)));
                try
                index = b(1);
                catch
                     disp a
                end
                positions = find(diff([1; freeSpace(:,i); 1])~=0)-1;
                lowerPoly(i) = positions(index);
                upperPoly(i) =  positions(index +1);
                
            end
            
            % divide freeSpace into numPoints blocks and find maximum
            
            blocks =  round(linspace(1,length(upperPoly),numPoints+1));
            blocks2 =  round(linspace(1,length(lowerPoly),numPoints+1));
            for i=1:numPoints
                [x,y] = min(upperPoly(blocks(i):blocks(i+1)));
                upperPoints(:,i) = [blocks(i)+y-1;x];
                [x,y] = max(lowerPoly(blocks2(i):blocks2(i+1)));
                lowerPoints(:,i) = [blocks2(i)+y-1;x];
            end
            
%             for i=1:numPoints
%                 [x,y] = min(upperPoly(blocks(i):blocks(i)));
%                 upperPoints(:,i) = [blocks(i)+y-1;x];
%                 [x,y] = max(lowerPoly(blocks2(i):blocks2(i)));
%                 lowerPoints(:,i) = [blocks2(i)+y-1;x];
%             end
            
                upperPoints = [(start-1+ upperPoints(1,:))*20;upperPoints(2,:) *20]';
                lowerPoints = [(start-1+ lowerPoints(1,:))*20;lowerPoints(2,:) *20]';

            
            
            
%             [~,upperPoly_binary] = max(freeSpace==0,[],1);
%             [~,lowerPoly_binary] = max(flip(freeSpace)==0,[],1);
            
            
            % fill gaps in polyline until reach numPoints
%             y1 = upperPoly_binary(1);
%             upperPoly(:,1) = [1;y1];
%             for i=2:length(upperPoly_binary)
%                 y2 = upperPoly_binary(i);
%                 if y2 ~= y1
%                     upperPoly = [upperPoly [i-1;y1]];
%                     upperPoly = [upperPoly [i;y2]];
%                 end
%                 y1 = y2;                
%             end
%             upperPoly = [upperPoly [length(upperPoly_binary);upperPoly_binary(end)]];
%             upperPoly
%             
%             isDone = 0;
%             gap = 1;
%             while ~isDone
%                 i = 0;
%                 while i <= length(upperPoly)-2
%                     
%                     if i == 0
%                         y1 = Inf;
%                     else
%                         y1 = upperPoly(2,i);
%                     end
%                     if i == length(upperPoly)-2
%                         y2 = Inf;
%                     else
%                         y2 = upperPoly(2,i+3);
%                     end
%                     
%                     if upperPoly(1,i+2)-upperPoly(1,i+1) <= gap && upperPoly(2,i+1) < y1  && upperPoly(2,i+1) < y2
%                         if i ~= 0
%                         upperPoly(1,i) = upperPoly(1,i)+1;
%                         end
%                         if i+3 <= length(upperPoly)
%                             upperPoly(1,i+3) = upperPoly(1,i+3)-1;
%                         end
%                         upperPoly(:,i+1:i+2) = [];
%                         if length(upperPoly) <= numPoints
%                             isDone = 1;
%                             break;
%                         end
%                         if i~=0 && i+1 <= length(upperPoly)
%                             if upperPoly(2,i) == upperPoly(2,i+1)
%                                 upperPoly(:,i:i+1) = [];
%                             end
%                         end
%                         if length(upperPoly) <= numPoints
%                             isDone = 1;
%                             break;
%                         end
%                         
%                         
%                     else
%                     i = i+2;
%                     end
% 
%                 end
%                 gap = gap + 1;
%             end
            
            
            
%             lowerPoly = 26-lowerPoly+1;
            
            
    
            
        end
    end
    
end

