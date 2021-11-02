classdef Ships
    
    properties
        
        numShips
        
        % constant
        id
        length
        width
        mass
        dir        
        maxPower = 1E6;
        
        % dynamic
        x
        x_utm
        y
        y_utm
        vx
        vy
        ax
        ay
        width_eff
        P
        ssq
        
        heading
        
        % objects
        latConPol
        lonConPol
        
        
        overtakeLevel
    end
    
    methods
        function o = Ships(river)
            numShips = 5;
            o.numShips = numShips;
            o.id = (1:numShips)';
            o.length = ones(numShips,1) * 100;
            o.width = ones(numShips,1) * 10;
            o.mass = ones(numShips,1) * 100E3;
            
             
            
            
            % init dynamic variables
            
%             o.x = 9000+[-2200 -2700 -4000 3700 4200 5500]'; %

            % Ship-following scenario
%             o.x = 9000-linspace(0,numShips*200, numShips)'; %
%             o.overtakeLevel = zeros(numShips,1);
%             o.dir = ones(numShips,1) ;
%             
            % 3 vs 3 scenario
%             o.x = 21000+[-2200 -2700 -3900 3600 4100 5300]'; %
%              o.overtakeLevel = [0 0 1 0 0 1]';
%              o.dir = [1 1 1 -1 -1 -1]' ;
             
             % 3 vs 2 scenario
%             o.x = 21000+[-2200 -2700 -3900 3600 4100 ]'; %
%              o.overtakeLevel = [0 0 1 0 0]';
%              o.dir = [1 1 1 -1 -1]' ;
            
            
%             o.x = 15500+[-2000 -2500 -3650 3500 4000 5150]'; %e Stelle im Rhein

            o.x = 40000-linspace(0,numShips*300, numShips)'; %
            o.overtakeLevel = [0 0 0 0 0]';
            o.dir = ones(numShips,1)*1;

            o.y = 150*ones(numShips,1) - 0*o.dir;
%             o.y(2) = -30;
%             o.y(6) = 150;
%            
            

            
            
            o.vx = 3*o.dir ;
            o.vy = zeros(numShips,1);
            
            o.ax = zeros(numShips,1);
            o.ay = zeros(numShips,1);
            o.width_eff = o.width;
            o.P = zeros(numShips,1);
            o.ssq = ones(numShips,1);
            o.heading.cf = zeros(numShips,1);
            o.heading.angle = zeros(numShips,1);
            o.heading.threePointsForRadiusCalc = cell(numShips,1);
            o.heading.box = cell(numShips,1);
            for i=1:numShips
                 o.heading.box{i} = o.createBox(i); 
            end
            o.x_utm = zeros(numShips,1);
            o.y_utm = zeros(numShips,1);
            
            
            o.latConPol = [];
            o.lonConPol = [];
            for i=1:numShips
                o.latConPol = [o.latConPol LatConPol(i)];
                o.lonConPol = [o.lonConPol LonConPol(i,o,river)];
            end
            
        end
        
        function o = simulateTimeStep(o,ID,SimSet,river)
            
            %%%%%%%% observations %%%%%%%%
            o.latConPol(ID) =  o.latConPol(ID).computeObs(o,river);
            o.lonConPol(ID) =  o.lonConPol(ID).computeObs(o,river);
            
            
            %%%%%%%% longitudinal %%%%%%%%
            if ID == 1 
                if mod(round(SimSet.t/SimSet.dT),200) == 1
                    o.P(ID) = (rand() * o.maxPower );
                    
                    if rand()< 0.1
                        o.P(ID) = 10000;
                    end
                end
                o.P(ID) = 3E5;
            elseif ID == 4
                o.P(ID) = 1E6;
            else
                Action = evaluatePolicy_lonConPol(o.lonConPol(ID).obs);
                o.P(ID) = max(0,o.maxPower * Action);
            end
%              o.P(ID) = 2E5;

             % o.P(ID) = 1E6; % set power for all vessels
             
             
             [o.lonConPol(ID),newAcc, newSquat, newCF] = o.lonConPol(ID).computeAccByPower(o,ID,SimSet,river);
            o.ax(ID) = o.dir(ID) * newAcc;
            o.ssq(ID) = newSquat;
            o.heading.cf(ID) = newCF;
           
            
            newVx= o.vx(ID) + o.ax(ID)*SimSet.dT;
            
%             if ID==1 && newVx < o.lonConPol(ID).obs_vStream + 7/3.6 
% %             if newVx < o.lonConPol(ID).obs_vStream + 7/3.6 
%                 newVx = o.lonConPol(ID).obs_vStream + 7/3.6 ;
%             end
            newX = o.x(ID) + 0.5 * (o.vx(ID) + newVx)*SimSet.dT;
            o.vx(ID) = newVx;
            o.x(ID) = newX;
            
            
            %%%%%%%% lateral %%%%%%%%
            Action = evaluatePolicy_latConPol(o.latConPol(ID).obs);
            
%             Action = Action /2;
            
            o.ay(ID) = Action * o.latConPol(ID).upperAccBound;
            newVy = o.vy(ID) + o.ay(ID) * SimSet.dT;
            if abs(newVy) > o.latConPol(ID).upperSpeedBound
                newVy = sign(newVy) * o.latConPol(ID).upperSpeedBound;
            end
            
            %  limit vy to max 0.1 vx !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if abs(newVy) > abs(o.vx(ID)/10)
                newVy = sign(newVy) *  abs(o.vx(ID)/10);
            end
            o.y(ID) = o.y(ID) + 0.5 * (o.vy(ID) + newVy)*SimSet.dT;
            o.vy(ID) = newVy;
            
            % update utm ship coordinates
            newPoint_utm = river.getPosUTMfromPos([o.x(ID);o.y(ID)]);
            o.x_utm(ID) = newPoint_utm(1);
            o.y_utm(ID) = newPoint_utm(2);
             o = o.computeHeadingFromCF(ID);
            
        end
        
        function [index, minDist] = findNextShipAhead(o,ID,num,overtakeLevel)
            temp = o.dir(ID) * (o.x + o.dir(ID)*o.length/2 - (o.x(ID)- o.dir(ID)*o.length(ID)/2)); % distance with lengths
            dist = o.dir(ID) * (o.x - o.x(ID));
            dist(temp <= 0 | o.dir ~= o.dir(ID) | o.id == ID |  ~ismember(o.overtakeLevel,overtakeLevel)) = inf;
            temp = sort(dist);
            minDist = temp(1:num);
            for i=1:num
                
                if minDist(i) == inf
                    index(i) = 0;
                else
                    index(i) = find(dist == minDist(i),1);
                end
            end
        end
        
        function [index, minDist] = findNextShipBehind(o,ID,num,overtakeLevel)
            temp = o.dir(ID) * ( o.x(ID) + o.dir(ID)*o.length(ID)/2 - (o.x - o.dir(ID)*o.length/2 ));
            dist = o.dir(ID) * ( o.x(ID) - o.x);
            dist(temp <= 0 | o.dir ~= o.dir(ID) | o.id == ID |  ~ismember(o.overtakeLevel,overtakeLevel)) = inf;
            temp = sort(dist);
            minDist = temp(1:num);
            for i=1:num
                
                if minDist(i) == inf
                    index(i) = 0;
                else
                    index(i) = find(dist == minDist(i),1);
                end
            end
        end
        
        function [index, minDist] = findOncomingShip(o,ID,num,overtakeLevel)
            temp = o.dir(ID) * (o.x + o.dir(ID).*o.length/2 - (o.x(ID)- o.dir(ID)*o.length(ID)/2)); % distance with lengths
            dist = o.dir(ID) * ( o.x - o.x(ID));
            dist(temp <= 0 | o.dir == o.dir(ID) | o.id == ID |  ~ismember(o.overtakeLevel,overtakeLevel)) = inf;
            temp = sort(dist);
            minDist = temp(1:num);
            for i=1:num
                
                if minDist(i) == inf
                    index(i) = 0;
                else
                    index(i) = find(dist == minDist(i),1);
                end
            end
            
        end
        
        function o = computeHeadingFromCF(o,ID)
            
            lastPoints = o.heading.threePointsForRadiusCalc{ID};
            lastPoints = [lastPoints [o.x_utm(ID);o.y_utm(ID)]];
            if size(lastPoints,2) ~=4
                headingAngle = 0;
            else
                lastPoints(:,4) = [];
                [Radius,centerPoint] = fit_circle_through_3_points( lastPoints');
                headingAngle = asin(o.length(ID) * (o.heading.cf(ID) - 0.5 )/ Radius);
                
                 if (centerPoint(1)-lastPoints(1,1))*(lastPoints(2,2)-lastPoints(2,1))-(centerPoint(2)-lastPoints(2,1)*(lastPoints(1,2)-lastPoints(1,1)))>0 % https://math.stackexchange.com/questions/274712/calculate-on-which-side-of-a-straight-line-is-a-given-point-located#:~:text=To%20see%20whether%20points%20on,point%20you%20are%20interested%20in.
                     headingAngle = - headingAngle;
                 end
                                
            end
%              headingAngle = 0; %%%%%%%%% REMOVE
            o.heading.angle(ID) = o.dir(ID) * headingAngle;
            o.heading.box{ID} = o.createBox(ID);
            o.heading.threePointsForRadiusCalc{ID} = lastPoints;
            
            
        end
        
        function box = createBox(o,ID)
             l = o.length(ID);
             w = o.width(ID);
           box = polyshape([-l/2, -l/2, l/2, l/2],[-w/2, w/2, w/2, -w/2]);
           o.heading.angle(ID) = (o.heading.angle(ID)+atan(o.vy(ID)/o.vx(ID)))/(2*pi)*360;
           box = rotate(box,o.heading.angle(ID));
%            box = rotate(box,atan(o.vy(ID)/o.vx(ID))/(2*pi)*360);
           box = translate(box,[o.x(ID) o.y(ID)]);
        
        end
        
        
        
    end
end

