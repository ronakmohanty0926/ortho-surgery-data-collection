function AnalyzeData(hapData, hapTime, forceData, forceTime)
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

tmp = 0;
for i = 1:length(sync_states)
    tmp = tmp + 1;
    drill_states(tmp,1) = sync_states(i).drill_states;
%     disp(drill_states(tmp,1));
end

% Compute All Parameters - Drill Tip Coordinates, Linear Speed, Angular
% Velocity, and Net Force
[drillt_coordinates, coordinates_series, sync_linearSpeeds, sync_netForces, timeHap, timeForce] = ComputeAllParam(device, sync_states, sync_forces);

% Segment Data
labelled_data = SegmentData(drill_states, sync_netForces, sync_linearSpeeds, drillt_coordinates, timeHap, timeForce);

len_approach = sum(~cellfun(@isempty,{labelled_data.position_approach}));
len_fapproach = sum(~cellfun(@isempty,{labelled_data.force_approach}));
len_rest = sum(~cellfun(@isempty,{labelled_data.position_rest}));
len_frest = sum(~cellfun(@isempty,{labelled_data.force_rest}));
len_cort1 = sum(~cellfun(@isempty,{labelled_data.position_cort1}));
len_fcort1 = sum(~cellfun(@isempty,{labelled_data.force_cort1}));
len_cort2 = sum(~cellfun(@isempty,{labelled_data.position_cort2}));
len_fcort2 = sum(~cellfun(@isempty,{labelled_data.force_cort2}));
len_retract = sum(~cellfun(@isempty,{labelled_data.position_retract}));
len_fretract = sum(~cellfun(@isempty,{labelled_data.force_retract}));

%Import struct data to local variables
count = 0;
while count < len_approach
    count = count + 1;
    position_approach(count,1) = labelled_data(count).position_approach;
    speed_approach(count,1) = labelled_data(count).speed_approach;
    timeHap_approach(count,1) = labelled_data(count).timeHap_approach;
end

count = 0;
while count < len_fapproach
    count = count + 1;
    force_approach(count,1) = labelled_data(count).force_approach;
    timeForce_approach(count,1) = labelled_data(count).timeForce_approach;
end

count = 0;
while count < len_rest
    count = count + 1;
    position_rest(count,1) = labelled_data(count).position_rest;
    speed_rest(count,1) = labelled_data(count).speed_rest;
    timeHap_rest(count,1) = labelled_data(count).timeHap_rest;
end

count = 0;
while count < len_frest
    count = count + 1;
    force_rest(count,1) = labelled_data(count).force_rest;
    timeForce_rest(count,1) = labelled_data(count).timeForce_rest;
end

count = 0;
while count < len_cort1
    count = count + 1;
    position_cort1(count,1) = labelled_data(count).position_cort1;
    speed_cort1(count,1) = labelled_data(count).speed_cort1;
    timeHap_cort1(count,1) = labelled_data(count).timeHap_cort1;
end

count = 0;
while count < len_fcort1
    count = count + 1;
    force_cort1(count,1) = labelled_data(count).force_cort1;
    timeForce_cort1(count,1) = labelled_data(count).timeForce_cort1;
end

count = 0;
while count < len_cort2
    count = count + 1;
    position_cort2(count,1) = labelled_data(count).position_cort2;
    speed_cort2(count,1) = labelled_data(count).speed_cort2;
    timeHap_cort2(count,1) = labelled_data(count).timeHap_cort2;
end    

count = 0;
while count < len_fcort2
    count = count + 1;
    force_cort2(count,1) = labelled_data(count).force_cort2;
    timeForce_cort2(count,1) = labelled_data(count).timeForce_cort2;
end

count = 0;
while count < len_retract
    count = count + 1;
    position_retract(count,1) = labelled_data(count).position_retract;
    speed_retract(count,1) = labelled_data(count).speed_retract;
    timeHap_retract(count,1) = labelled_data(count).timeHap_retract;
end    

count = 0;
while count < len_fretract
    count = count + 1;
    force_retract(count,1) = labelled_data(count).force_retract;
    timeForce_retract(count,1) = labelled_data(count).timeForce_retract;
end

f1 = figure('units','normalized','outerposition',[0 0 1 1]);
xlabel('Time in Seconds');
ylabel('Force in Newtons');
hold on
% plot(force_approach(1:end,2), force_approach(1:end,1),'LineWidth',2,'Color','k');
% plot(force_rest(1:end,2), force_rest(1:end,1),'LineWidth',2,'Color','k');
plot(timeForce_cort1(1:end), force_cort1(1:end),'LineWidth',2,'Color','r');
plot(timeForce_cort2(1:end), force_cort2(1:end),'LineWidth',2,'Color','g');
plot(timeForce_retract(1:end), force_retract(1:end),'LineWidth',2,'Color','b');
hold off
legend('First Cortex', 'Second Cortex','Retract');  
% f1 = plot(timeForce(1:end), sync_netForces(1:end),'LineWidth',2,'Color','r');
% xlim([force_cort1(1,2) force_retract(end,2)]);
% ylim([force_cort1(1,1) force_retract(end,1)]);
% xticks(force_cort1(1,2):2:force_retract(end,2));
% yticks(force_cort1(1,1):5:force_retract(end,1));
saveas(f1,'HB4_Force.png','png');

f2 = figure('units','normalized','outerposition',[0 0 1 1]);
xlabel('Time in Seconds');
ylabel('Speed in mm/s');
hold on
% plot(speed_approach(1:end,2), speed_approach(1:end,1),'LineWidth',2,'Color','k');
% plot(speed_rest(1:end,2), speed_rest(1:end,1),'LineWidth',2,'Color','k');
plot(timeHap_cort1(1:end), speed_cort1(1:end),'LineWidth',2,'Color','r');
plot(timeHap_cort2(1:end), speed_cort2(1:end),'LineWidth',2,'Color','g');
plot(timeHap_retract(1:end), speed_retract(1:end),'LineWidth',2,'Color','b');
hold off
legend('First Cortex', 'Second Cortex','Retract');  
% xlim([speed_cort1(1,2) speed_retract(end,2)]);
% ylim([speed_cort1(1,1) speed_retract(end,1)]);
% xticks(speed_cort1(1,2):2:speed_retract(end,2));
% yticks(speed_cort1(1,1):5:speed_retract(end,1));
saveas(f2,'HB4_Speed.png','png');

f3 = figure('units','normalized','outerposition',[0 0 1 1]);
xlabel('Time in Seconds');
ylabel('Stylus Tip Position in mm');
hold on
% plot(position_approach(1:end,2), position_approach(1:end,1),'LineWidth',2,'Color','k');
% plot(position_rest(1:end,2), position_rest(1:end,1),'LineWidth',2,'Color','k');
plot(timeHap_cort1(1:end), position_cort1(1:end),'LineWidth',2,'Color','r');
plot(timeHap_cort2(1:end), position_cort2(1:end),'LineWidth',2,'Color','g');
plot(timeHap_retract(1:end), position_retract(1:end),'LineWidth',2,'Color','b');
hold off
legend('First Cortex', 'Second Cortex','Retract');       
% xlim([position_cort1(1,2) position_retract(end,2)]);
% xticks(position_cort1(1,1):2:position_retract(end,1));
saveas(f3,'HB4_Position.png','png');
close('all');


% SaveAllParamPlotsToGif(device, sync_states, sync_forces, 'YBUser07_0.gif');
% SaveDeviceStatesToGif(device, sync_states, 'fancy', 50, '65Trial_Drill1.gif');

% clear;