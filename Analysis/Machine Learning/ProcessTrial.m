function [drillPosition, sync_linearSpeeds, sync_acceleration, sync_jerk, sync_speedRatio, netForces, timeHap, newtimeForce, drill_states] = ProcessTrial(hapData, hapTime, forceData, forceTime)
% function AnalyzeData(hapData, hapTime, forceData, forceTime, tempData)
% clear;
% Read Log Files
% hapticsDataTable = load('Log_User06_2.txt');
% forceDataTable = load('FTUser06_2.txt');
% temperatureDataTable = fopen('TTUser06_2.txt');
hapticsDataTable = load(hapData);
forceDataTable = load(forceData);

% Read Time Stamp files
hapticsTimeStamp = fopen(hapTime);
forceTimeStamp = fopen(forceTime);

% Read Log Data into Struct
[device, hapticsStates] = ReadHapticsData(hapticsDataTable);
forces = ReadForceTorqueData(forceDataTable);
% [temperature, thermal_times] = ReadTemperatureData(temperatureDataTable);

% Read Time Stamp Data into Struct
hapticsTimes = ReadHapticsTimeStamp(hapticsTimeStamp);
forceTimes = ReadForceTimeStamp(forceTimeStamp);

% Synchronize Data
[sync_states, sync_forces] = SyncData(hapticsTimes, hapticsStates, forceTimes, forces);

count = 0;
for i = 1:length(sync_states)
    count = count + 1;
    drill_states(count,1) = sync_states(i).drill_states;
end

% Compute All Parameters - Drill Tip Coordinates, Linear Speed, Angular
% Velocity, and Net Force
[drillPosition, sync_linearSpeeds, sync_acceleration, sync_jerk, sync_speedRatio, sync_netForces, timeHap, timeForce] = ComputeAllParam(device, sync_states, sync_forces);


% Sync data instances from Force with Haptics data
% count = 1;
% skipFrame = 3;
% idx = 1;
% netForces(1,1) = sync_netForces(1,1);
% newtimeForce(1,1) = timeForce(1,1);

idx = 1;
  
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
% while count < (length(sync_netForces) - skipFrame)
%     count = count + skipFrame;
%     idx = idx + 1;
%     netForces(idx,1) = sync_netForces(count,1);
%     newtimeForce(idx,1) = timeForce(count,1);
% end

%Compute Delta for each Data
% for i = 2:length(drillt_coordinates)
%     delta_position(i,1) = drillt_coordinates(i,1) - drillt_coordinates(i-1,1);
%     delta_speed(i,1) = sync_linearSpeeds(i,1) - sync_linearSpeeds(i-1,1);
%     delta_force(i,1) = netForces(i,1) - netForces(i-1,1);
% end
% delta_position(1,1) = 0;
% delta_speed(1,1) = 0;
% delta_force(1,1) = 0;

