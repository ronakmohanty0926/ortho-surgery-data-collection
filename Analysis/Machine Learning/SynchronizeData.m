function [sync_states, sync_forces, sync_torques, sync_temperatures] = SynchronizeData(hapticsTimes, hapticsStates, forceTimes, forces, torques, temperature, thermal_times)
% SynchronizeData synchronizes the OpenHaptics API raw data and F/T sensor
% data. The assumption here is that the OpenHaptics device recording is 
% started and ended first. Therefore, there is a time lag between the start 
% and end of the F/T data log with respect the OpenHaptics data.
% Load the haptics and forces raw data and time structs into SynchronizeData 
% to synchronize these two data logs
%
% e.g. [sync_states, sync_forces, sync_torques] = SynchronizeData(hapticsTimes, hapticsStates, forceTimes, forces, torques)

hapticsDataLen = length(hapticsStates);
forceDataLen = length(forces);
thermalDataLen = length(temperature);

%Start Time Sync for Haptics and Thermal Camera

if hapticsTimes(1).start_second > thermal_times(1).start_second
    thermal_times(1).start_minute = thermal_times(1).start_minute - 1;
    thermal_times(1).start_second = thermal_times(1).start_second + 60;
end
hap_diff_second = thermal_times(1).start_second - hapticsTimes(1).start_second;

if hapticsTimes(1).start_minute > thermal_times(1).start_minute
    thermal_times(1).start_hour = thermal_times(1).start_hour - 1;
    thermal_times(1).start_minute = thermal_times(1).start_minute + 60;
end
hap_diff_minute = thermal_times(1).start_minute - hapticsTimes(1).start_minute;
hap_diff_hour = thermal_times(1).start_hour - hapticsTimes(1).start_hour; 

haptics_lag_duration = (hap_diff_hour*3600) + (hap_diff_minute*60) + hap_diff_second;

temp_time = 0;
hapframeCount = 0;

for i = 1: hapticsDataLen
    if temp_time < haptics_lag_duration
        hapframeCount = i;
        temp_time = hapticsStates(i).time_stamp;
    end
end
hapticsStartFrame = hapframeCount + 1;
frameLen = hapticsDataLen - hapframeCount;
% frameLen = 1000;

if hapticsTimes(1).end_second > thermal_times(1).end_second
    thermal_times(1).end_minute = thermal_times(1).end_minute - 1;
    thermal_times(1).end_second = thermal_times(1).end_second + 60;
end
thermal_end_diff_second = thermal_times(1).end_second - hapticsTimes(1).end_second;

if hapticsTimes(1).end_minute > thermal_times(1).end_minute
    thermal_times(1).end_hour = thermal_times(1).end_hour - 1;
    thermal_times(1).end_minute = thermal_times(1).end_minute + 60;
end
thermal_end_diff_minute = thermal_times(1).end_minute - hapticsTimes(1).end_minute;
thermal_end_diff_hour = thermal_times(1).end_hour - hapticsTimes(1).end_hour; 

thermal_end_lag_duration = (thermal_end_diff_hour*3600) + (thermal_end_diff_minute*60) + thermal_end_diff_second;

thDataLen = thermalDataLen - (8.*thermal_end_lag_duration);

%Start Time Sync for Force and Thermal Camera

if forceTimes(1).start_second > thermal_times(1).start_second
    thermal_times(1).start_minute = thermal_times(1).start_minute - 1;
    thermal_times(1).start_second = thermal_times(1).start_second + 60;
end
force_diff_second = thermal_times(1).start_second - forceTimes(1).start_second;

if forceTimes(1).start_minute > thermal_times(1).start_minute
    thermal_times(1).start_hour = thermal_times(1).start_hour - 1;
    thermal_times(1).start_minute = thermal_times(1).start_minute + 60;
end
force_diff_minute = thermal_times(1).start_minute - forceTimes(1).start_minute;
force_diff_hour = thermal_times(1).start_hour - forceTimes(1).start_hour; 

force_lag_duration = (force_diff_hour*3600) + (force_diff_minute*60) + force_diff_second;

forceStartFrame = (1000.*force_lag_duration) + 1;

if hapticsTimes(1).end_second > forceTimes(1).end_second
    forceTimes(1).end_minute = forceTimes(1).end_minute - 1;
    forceTimes(1).end_second = forceTimes(1).end_second + 60;
end
force_end_diff_second = forceTimes(1).end_second - hapticsTimes(1).end_second;

if hapticsTimes(1).end_minute > forceTimes(1).end_minute
    forceTimes(1).end_hour = forceTimes(1).end_hour - 1;
    forceTimes(1).end_minute = forceTimes(1).end_minute + 60;
end
force_end_diff_minute = forceTimes(1).end_minute - hapticsTimes(1).end_minute;
force_end_diff_hour = forceTimes(1).end_hour - hapticsTimes(1).end_hour; 

force_end_lag_duration = (force_end_diff_hour*3600) + (force_end_diff_minute*60) + force_end_diff_second;
forceSkipFrames = 1000.*force_end_lag_duration;

forceStartFrame = (1000.*force_lag_duration) + 1;
fDataLen = forceDataLen - (1000.*force_end_lag_duration);
forceIntervalTime = 0.001;
force_count = 0;
force_time = 0;

for j = forceStartFrame:fDataLen
    force_count = force_count + 1;
    force_time = force_time + forceIntervalTime;
    sync_forces(force_count).x_axis = forces(j).x_axis;
    sync_forces(force_count).y_axis = forces(j).y_axis;
    sync_forces(force_count).z_axis = forces(j).z_axis;
    sync_torques(force_count).x_axis = torques(j).x_axis;
    sync_torques(force_count).y_axis = torques(j).y_axis;
    sync_torques(force_count).z_axis = torques(j).z_axis;
    sync_forces(force_count).time_diff = force_time;    
end

thermalIntervalTime = 1./8;
thermal_time = 0;

for k = 1:thDataLen
    thermal_time = thermal_time + thermalIntervalTime;
    sync_temperatures(k).cortex1 = temperature(k).cortex1;
    sync_temperatures(k).cortex2 = temperature(k).cortex2;
    sync_temperatures(k).time_diff = thermal_time;
end

hap_count = 0;

for m = (hapticsStartFrame):hapticsDataLen
    hap_count = hap_count + 1;
    sync_states(hap_count).base_angle = hapticsStates(m).base_angle;
    sync_states(hap_count).angle_link01 = hapticsStates(m).angle_link01;
    sync_states(hap_count).angle_link12 = hapticsStates(m).angle_link12;
    sync_states(hap_count).stylus_gamma = hapticsStates(m).stylus_gamma;
    sync_states(hap_count).stylus_beta = hapticsStates(m).stylus_beta;
    sync_states(hap_count).stylus_alpha = hapticsStates(m).stylus_alpha;    
    sync_states(hap_count).time_stamp = hapticsStates(m).time_stamp;
end

for n = 1:length(sync_states)
    sync_states(n).time_diff = sync_states(n).time_stamp - sync_states(1).time_stamp;
end


    





