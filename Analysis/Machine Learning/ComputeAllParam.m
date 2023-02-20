function [drillPosition, sync_linearSpeeds, sync_acceleration, sync_jerk, sync_speedRatio, sync_netForces, timeHap, timeForce] = ComputeAllParam(device, sync_states, sync_forces)


coordinates_series = ComputeCoordinatesForStates(device, sync_states);

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
sync_acceleration = ComputeAcceleration(sync_linearSpeeds, timeHap);
sync_jerk = ComputeJerk(sync_acceleration, timeHap);
sync_speedRatio = ComputeSpeedRatio(sync_linearSpeeds);
sync_netForces = ComputeNetForces(sync_forces);
count = 0;

for j = 1:length(sync_forces)
    timeForce(j,1) = sync_forces(j).time_diff;
end