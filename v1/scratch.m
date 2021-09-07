utm = [];
for i=1:5
    box = o.heading.box{i}.Vertices;
    for j=1:4        
        utm = [utm; river.getPosUTMfromPos([box(j,1) box(j,2)])']; 
    end
end

for i=1:7
    utm = [utm;river.getPosUTMfromPos(o.latConPol(3).basePointsL(i,:))'];
end

for i=1:7
    utm = [utm;river.getPosUTMfromPos(o.latConPol(3).basePointsU(i,:))'];
end