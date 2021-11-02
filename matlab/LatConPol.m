classdef LatConPol
    
    properties
        
        % constant
        id
        safetyDist
        upperSpeedBound
        upperTTCbound
        lowerTTCbound
        upperRiverBound
        upperAccBound
        numFeatures
        lookAheadDistanceFreeSpaceArea
        test
        
        %dynamic
        obs
        basePointsL
        basePointsU
        
        obs_ttcbasePointsL
        obs_ybasePointsL
        obs_ttcbasePointsU
        obs_ybasePointsU
        TTCvectorUp
        TTCvectorLow
        POSvectorUp
        POSvectorLow
        
        
        obs_ttcShipAhead
        obs_yShipAhead
        
        obs_ttcShipBehind
        obs_yShipBehind
        
        obs_ttcOncomingShip
        obs_yOncomingShip
    end
    
    methods
        function o = LatConPol(id)
            o.id = id;
            o.upperSpeedBound = 4;
            o.upperTTCbound = 1000;
            o.lowerTTCbound = -1000;
            o.upperRiverBound = 500;
            o.upperAccBound = 0.01;
            o.safetyDist = 0;
            o.numFeatures = 14;
            o.lookAheadDistanceFreeSpaceArea = 4500;
            
        end
        
        function o = computeObs(o,ships, river)
            ID = o.id;
            dir = ships.dir(ID);
            num = 7;
            
            
            % river basepoints
            
            % real observation
            [ o.basePointsL ,  o.basePointsU ] = river.detectFreeSpace(ID,ships,o.lookAheadDistanceFreeSpaceArea, num);
            %o.basePointsU = o.basePointsU - 0;
            %o.basePointsL = o.basePointsU - 40;
            o.obs_ttcbasePointsL = dir*[(o.basePointsL(:,1) - ships.x(ID) )/abs(ships.vx(ID))...
                (o.basePointsL(:,1) - ships.x(ID) )/abs(ships.vx(ID))];
            o.obs_ybasePointsL = [o.basePointsL(:,2) - ships.y(ID) o.basePointsL(:,2) - ships.y(ID)];
            
            o.obs_ttcbasePointsU = dir*[(o.basePointsU(:,1) - ships.x(ID) )/abs(ships.vx(ID))...
                (o.basePointsU(:,1) - ships.x(ID) )/abs(ships.vx(ID))];
            o.obs_ybasePointsU = [o.basePointsU(:,2) - ships.y(ID) o.basePointsU(:,2) - ships.y(ID)];
            
            
           % artificial coastline
%             o.basePointsL = river.getClosestBasePointsL(ships.x(ID),dir,num)';
%             o.basePointsU = river.getClosestBasePointsU(ships.x(ID),dir,num)';
%             
%             o.obs_ttcbasePointsL = dir*[(o.basePointsL(:,1) - ships.x(ID) - ships.length(ID)/2)/abs(ships.vx(ID))...
%                 (o.basePointsL(:,1) - ships.x(ID) + ships.length(ID)/2)/abs(ships.vx(ID))];
%             o.obs_ybasePointsL = [o.basePointsL(:,2) - ships.y(ID) o.basePointsL(:,2) - ships.y(ID)];
%             
%             o.obs_ttcbasePointsU = dir*[(o.basePointsU(:,1) - ships.x(ID) - ships.length(ID)/2)/abs(ships.vx(ID))...
%                 (o.basePointsU(:,1) - ships.x(ID) + ships.length(ID)/2)/abs(ships.vx(ID))];
%             o.obs_ybasePointsU = [o.basePointsU(:,2) - ships.y(ID) o.basePointsU(:,2) - ships.y(ID)];
%             
%             
            % Ship ahead observations
            num = 3;
            [shipAheadID, distShipAhead] = ships.findNextShipAhead(ID,num,(0:ships.overtakeLevel(ID)-1));
            for i=1:num
                if shipAheadID(i)  ~= 0
                    if abs(ships.vx(ID)) > abs(ships.vx(shipAheadID(i)))
                        if dir == 1
                            o.obs_ttcShipAhead(i,:) = [min(ships.heading.box{shipAheadID(i)}.Vertices(1:2,1)) - max(ships.heading.box{ID}.Vertices(3:4,1)), max(ships.heading.box{shipAheadID(i)}.Vertices(3:4,1)) - min(ships.heading.box{ID}.Vertices(1:2,1))]/abs(ships.vx(ID) - ships.vx(shipAheadID(i)));
                            o.obs_yShipAhead(i,:) = [ships.heading.box{shipAheadID(i)}.Vertices(2,2) - ships.heading.box{ID}.Vertices(4,2), ships.heading.box{shipAheadID(i)}.Vertices(3,2) - ships.heading.box{ID}.Vertices(1,2)];
                        else
                            o.obs_ttcShipAhead(i,:) = [min(ships.heading.box{ID}.Vertices(1:2,1)) - max(ships.heading.box{shipAheadID(i)}.Vertices(3:4,1)), max(ships.heading.box{ID}.Vertices(3:4,1)) - min(ships.heading.box{shipAheadID(i)}.Vertices(1:2,1))]/abs(ships.vx(ID) - ships.vx(shipAheadID(i)));
                            o.obs_yShipAhead(i,:) = [ships.heading.box{shipAheadID(i)}.Vertices(4,2) - ships.heading.box{ID}.Vertices(2,2), ships.heading.box{shipAheadID(i)}.Vertices(1,2) - ships.heading.box{ID}.Vertices(3,2)];
                            
                        end
                    else
                        o.obs_ttcShipAhead(i,:) = [inf;inf];
                        o.obs_yShipAhead(i,:) = [0;0];
                    end
                else
                    o.obs_ttcShipAhead(i,:) =  [inf;inf];
                    o.obs_yShipAhead(i,:) = [0;0];
                end
            end
            
            % ship behind observations
            num = 3;
            [shipBehindID, distShipBehind] = ships.findNextShipBehind(ID,num,ships.overtakeLevel(ID)+1:3);
            for i=1:num
                if shipBehindID(i)  ~= 0
                    if abs(ships.vx(ID)) < abs(ships.vx(shipBehindID(i))) 
                         if dir == 1
                            o.obs_ttcShipBehind(i,:) = [min(ships.heading.box{ID}.Vertices(1:2,1)) - max(ships.heading.box{shipBehindID(i)}.Vertices(3:4,1)), max(ships.heading.box{ID}.Vertices(3:4,1)) - min(ships.heading.box{shipBehindID(i)}.Vertices(1:2,1))]/abs(ships.vx(ID) - ships.vx(shipBehindID(i)));
                            o.obs_yShipBehind(i,:) = [ships.heading.box{shipBehindID(i)}.Vertices(4,2) - ships.heading.box{ID}.Vertices(2,2), ships.heading.box{shipBehindID(i)}.Vertices(1,2) - ships.heading.box{ID}.Vertices(3,2) ];
                        else
                            o.obs_ttcShipBehind(i,:) = [min(ships.heading.box{shipBehindID(i)}.Vertices(1:2,1)) - max(ships.heading.box{ID}.Vertices(3:4,1)), max(ships.heading.box{shipBehindID(i)}.Vertices(3:4,1)) - min(ships.heading.box{ID}.Vertices(1:2,1))]/abs(ships.vx(ID) - ships.vx(shipBehindID(i)));
                            o.obs_yShipBehind(i,:) = [ships.heading.box{shipBehindID(i)}.Vertices(2,2) - ships.heading.box{ID}.Vertices(4,2), ships.heading.box{shipBehindID(i)}.Vertices(3,2) - ships.heading.box{ID}.Vertices(1,2) ];
                            
                        end
                    else
                        o.obs_ttcShipBehind(i,:) = [inf;inf] ;
                        o.obs_yShipBehind(i,:)   = [0;0];
                    end
                else
                    o.obs_ttcShipBehind(i,:) = [inf;inf];
                    o.obs_yShipBehind(i,:)  = [0;0];
                end
            end
            
            
            % oncoming ship observations
            num = 3;
            [oncomingShipID, oncomingShipDist] = ships.findOncomingShip(ID,num,0:3);
            for i=1:num
                
                if oncomingShipID(i)  ~= 0
                    if dir == 1
                        o.obs_ttcOncomingShip(i,:) = [min(ships.heading.box{oncomingShipID(i)}.Vertices(1:2,1)) - max(ships.heading.box{ID}.Vertices(3:4,1)), max(ships.heading.box{oncomingShipID(i)}.Vertices(3:4,1)) - min(ships.heading.box{ID}.Vertices(1:2,1))]/abs(ships.vx(ID) - ships.vx(oncomingShipID(i)));
                        o.obs_yOncomingShip(i,:) = [ships.heading.box{oncomingShipID(i)}.Vertices(1,2) - ships.heading.box{ID}.Vertices(3,2), ships.heading.box{oncomingShipID(i)}.Vertices(4,2) - ships.heading.box{ID}.Vertices(2,2) ];
                    else
                        o.obs_ttcOncomingShip(i,:) = [min(ships.heading.box{ID}.Vertices(1:2,1)) - max(ships.heading.box{oncomingShipID(i)}.Vertices(3:4,1)), max(ships.heading.box{ID}.Vertices(3:4,1)) - min(ships.heading.box{oncomingShipID(i)}.Vertices(1:2,1))]/abs(ships.vx(ID) - ships.vx(oncomingShipID(i)));
                        o.obs_yOncomingShip(i,:) = [ships.heading.box{oncomingShipID(i)}.Vertices(3,2) - ships.heading.box{ID}.Vertices(1,2), ships.heading.box{oncomingShipID(i)}.Vertices(2,2) - ships.heading.box{ID}.Vertices(4,2) ];
                        
                    end
%                     o.obs_ttcOncomingShip(i,:)= [(oncomingShipDist(i)-ships.length(oncomingShipID(i))/2-ships.length(ID)/2)/abs(ships.vx(oncomingShipID(i))-ships.vx(ID));
%                         (oncomingShipDist(i)+ships.length(oncomingShipID(i))/2+ships.length(ID)/2)/abs(ships.vx(oncomingShipID(i))-ships.vx(ID))];
%                     o.obs_yOncomingShip(i,:)  = [1;1]*ships.y(oncomingShipID(i)) - ships.y(ID);
%                     width_oncomingShip = ships.width(oncomingShipID(i));
%                     
                else
                    o.obs_ttcOncomingShip(i,:) = [inf;inf];
                    o.obs_yOncomingShip(i,:)  = [0;0];
                end
            end
            
            
            try
            if dir == 1
                [o.TTCvectorUp,o.POSvectorUp] = o.computeTTCvector('Up',o.numFeatures/2,ships,...
                    o.obs_ttcbasePointsU,...
                    [o.obs_ttcShipBehind;
                    o.obs_ttcOncomingShip],...
                    o.obs_ybasePointsU,...
                    [o.obs_yShipBehind;
                    o.obs_yOncomingShip]);
                
                [o.TTCvectorLow,o.POSvectorLow] = o.computeTTCvector('Low',o.numFeatures/2,ships,...
                    o.obs_ttcbasePointsL,...
                    o.obs_ttcShipAhead,...
                    o.obs_ybasePointsL,...
                    o.obs_yShipAhead );
                
                
            else
                [o.TTCvectorUp,o.POSvectorUp] = o.computeTTCvector('Up',o.numFeatures/2,ships,...
                    o.obs_ttcbasePointsU,...
                    o.obs_ttcShipAhead,...
                    o.obs_ybasePointsU,...
                    o.obs_yShipAhead);
                
                [o.TTCvectorLow,o.POSvectorLow] = o.computeTTCvector('Low',o.numFeatures/2,ships,...
                    o.obs_ttcbasePointsL,...
                    [o.obs_ttcShipBehind;
                    o.obs_ttcOncomingShip],...
                    o.obs_ybasePointsL,...
                    [o.obs_yShipBehind ;
                    o.obs_yOncomingShip]);
                
            end
            catch
                disp BUG!!!
            end
            
            % observation
            o.obs = [
                ships.ay(ID)/o.upperAccBound...
                ships.vy(ID)/o.upperSpeedBound...
                o.TTCvectorUp'/o.upperTTCbound...
                o.TTCvectorLow'/o.upperTTCbound...
                o.POSvectorUp'/(2*o.upperRiverBound)...
                o.POSvectorLow'/(2*o.upperRiverBound)
                ];
            
               o.obs = [
                ships.ay(ID)/o.upperAccBound...
                ships.vy(ID)/o.upperSpeedBound...
                o.TTCvectorUp'/o.upperTTCbound...
                o.TTCvectorLow'/o.upperTTCbound...
                tanh(o.POSvectorUp'/300)...
                tanh(o.POSvectorLow'/300)
                ];
        end
        
        
        
        
        function [TTCvector, POSvector] = computeTTCvector(o,direction,numEntries,ships,c_ttc, v_ttc, c_pos,  v_pos)
            
            pos =  [c_pos;  v_pos];
            ttc = [c_ttc; v_ttc];
            
            if strcmp(direction,'Up')
                sortDir = 'ascend';
            else
                sortDir = 'descend';
            end
            [~, sortInd] = sort(pos(:,1),sortDir);
            pos = pos(sortInd,:);
            ttc = ttc(sortInd,:);
            
            toDel = zeros(length(pos),2);
            for i=1:length(pos)-1
                ttcLow = ttc(i,1);
                ttcHigh = ttc(i,2);
                for j=i+1:length(pos)
                    if ttc(j,1) > ttcLow && ttc(j,2) < ttcHigh
                        % tuple is completely hidden
                        ttc(j,1) = ttcLow;
                        ttc(j,2) = ttcHigh;
                        toDel(j,1) = 1;
                        toDel(j,2) = 1;
                    elseif ttc(j,1)< ttcLow && ttc(j,2) < ttcHigh && ttc(j,2) > ttcLow
                        ttc(j,2) = ttcLow;
                        toDel(j,2) = 1;
                    elseif ttc(j,1) > ttcLow && ttc(j,2) > ttcHigh && ttc(j,1)  < ttcHigh
                        ttc(j,1) = ttcHigh;
                        toDel(j,1) = 1;
                    end
                end
            end
            
            ttc = [ttc(:,1);ttc(:,2)];
            pos = [pos(:,1);pos(:,2)];
            toDel = [toDel(:,1);toDel(:,2)];
            
            ttc(boolean(toDel)) = [];
            pos(boolean(toDel)) = [];
            
            [ttc, sortInd] = sort(ttc);
            pos = pos(sortInd);
            
            % delete double entries
             toDel = diff(ttc)==0 & diff(pos)==0;
             ttc(toDel) = [];
             pos(toDel) = [];
            
            % too many negative ttc entries -> delete
            while ttc(2)<0
                ttc(1) = [];
                pos(1) = [];
            end
            
            if length(ttc) < numEntries
                ttc(end+1:numEntries) = o.upperTTCbound;
                pos(end+1:numEntries) = 0;
            end          
            
            
            if strcmp(direction,'Up')
                pos = pos - o.safetyDist ;
            else
                pos = pos + o.safetyDist ;
            end
            TTCvector = clip(ttc(1:numEntries),o.lowerTTCbound, o.upperTTCbound);
            POSvector = clip(pos(1:numEntries),-o.upperRiverBound,o.upperRiverBound);
            
        end
    end
end

