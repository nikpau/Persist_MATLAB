% run computeVesselFairway.m first
tic
dataToPlot = 5; % 4: stream velocity in x
                % 5: stream velocity in y
                % 6: water depth
                % 7: abs stream velocity
                
%load inp file
opts=detectImportOptions('MW_grid_500m_655_665_UTM32.inp','FileType','text');
data2 = readmatrix('MW_grid_500m_655_665_UTM32.inp',opts);

index_ones = find(data2(:,1)==1);
PointCoordinates = data2(index_ones(1):index_ones(2)-9,2:3); %Matrix with all point coordinates
PointDataset = data2(index_ones(2):end,2:7); % Matrix with all point dataset (like stream velocity,...)



% 
% figure
%  scatter(PointCoordinates(:,1),PointCoordinates(:,2)) 
%  daspect([1 1 1])
%  xlabel('km')
%  ylabel('km')
%  
 

    streamVel = zeros(length(PointCoordinates),1);
 for i=1:length(PointCoordinates)
 streamVel(i) = sqrt(PointDataset(i,1)^2 + PointDataset(i,2)^2);
 end

 streamVel  = reshape(streamVel,26,[]);
 PointDataset = reshape(PointDataset,26,[],6);
 PointCoordinates = reshape(PointCoordinates,26,[],2);
 toc
 %save river data to file
%  writematrix(gridPoints2,'Flusskarte_gebogen_km642_km647.csv' );
