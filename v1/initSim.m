function [river, ships, SimSet, viz] = initSim()

SimSet.dT = 1;
SimSet.vizInUTM = 1;

SimSet.isDone = 0;
SimSet.t = 0;

river = River();
ships = Ships(river);



viz = Viz(ships.numShips,river,ships,SimSet);

end