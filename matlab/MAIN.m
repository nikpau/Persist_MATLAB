tic

[river, ships, SimSet, viz] = initSim();

while ~SimSet.isDone
    
    for shipID = 1:ships.numShips
        
        
        ships = ships.simulateTimeStep(shipID, SimSet, river);
        
    end
    SimSet.t = SimSet.t + SimSet.dT;
    [viz, SimSet] = viz.visualize(ships,river,SimSet);
    
    
     if SimSet.t == 900
%         SimSet.isDone = 1;
    end
end

toc