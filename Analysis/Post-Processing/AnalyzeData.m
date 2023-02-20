function [drillPosition, sync_linearSpeeds, sync_netForces, timeHap, timeForce, drill_states] = AnalyzeData(hapData, hapTime, forceData, forceTime)
% function AnalyzeData(hapData, hapTime, forceData, forceTime, tempData)
% clear;
% Read Log Files
% hapticsDataTable = load('Log_User06_2.txt');
% forceDataTable = load('FTUser06_2.txt');
% temperatureDataTable = fopen('TTUser06_2.txt');
hapticsDataTable = load(hapData);
forceDataTable = load(forceData);
% temperatureDataTable = fopen(forceData);

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
%     disp(drill_states(tmp,1));
end

% Compute All Parameters - Drill Tip Coordinates, Linear Speed, Angular
% Velocity, and Net Force
[drillPosition, sync_linearSpeeds, sync_netForces, timeHap, timeForce] = ComputeAllParam(device, sync_states, sync_forces);

