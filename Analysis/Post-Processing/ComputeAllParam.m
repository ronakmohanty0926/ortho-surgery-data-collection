function [drillPosition, sync_linearSpeeds, netForces, timeHap, newtimeForce] = ComputeAllParam(device, sync_states, sync_forces)

coordinates_series = ComputeCoordinatesForStates(device, sync_states);
sync_linearSpeeds = ComputeLinearSpeeds(device, sync_states);
sync_netForces = ComputeNetForces(sync_forces);

for j = 1:length(sync_states)   
    drillPosition(j,1) = coordinates_series(j).drilltip(1);
    drillPosition(j,2) = coordinates_series(j).drilltip(2);
    drillPosition(j,3) = coordinates_series(j).drilltip(3);
end

for i = 1:length(sync_states)
    timeHap(i,1) = sync_states(i).time_diff;
%     drillt_coordinates(i,1) = coordinates_series(i).drilltip(2);
end

sync_linearSpeeds = ComputeLinearSpeeds(drillPosition, timeHap);
sync_netForces = ComputeNetForces(sync_forces);
% count = 0;

for j = 1:length(sync_forces)
    timeForce(j,1) = sync_forces(j).time_diff;
end

idx = 0;
  
for j = 1:3:length(sync_netForces) 
  if j < (length(sync_netForces)-3)
      idx = idx + 1;
      netForces(idx,1) = sync_netForces(j,1);
      newtimeForce(idx,1) = timeForce(j,1);
      j_end = j;
  end
end

for k = j_end:length(sync_netForces) 
    idx = idx + 1;
    netForces(idx,1) = sync_netForces(k,1);
    newtimeForce(idx,1) = timeForce(k,1);
end

