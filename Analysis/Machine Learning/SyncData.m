function [sync_states, sync_forces] = SyncData(hapticsTimes, hapticsStates, forceTimes, forces)
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

%Start Time Sync for Haptics and Thermal Camera

if hapticsTimes(1).start_second > forceTimes(1).start_second
    forceTimes(1).start_minute = forceTimes(1).start_minute - 1;
    forceTimes(1).start_second = forceTimes(1).start_second + 60;
end
hap_diff_second = forceTimes(1).start_second - hapticsTimes(1).start_second;

if hapticsTimes(1).start_minute > forceTimes(1).start_minute
    forceTimes(1).start_hour = forceTimes(1).start_hour - 1;
    forceTimes(1).start_minute = forceTimes(1).start_minute + 60;
end
hap_diff_minute = forceTimes(1).start_minute - hapticsTimes(1).start_minute;
hap_diff_hour = forceTimes(1).start_hour - hapticsTimes(1).start_hour; 

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
% frameLen = hapticsDataLen - hapframeCount;

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

fDataLen = forceDataLen - forceSkipFrames;
force_count = 0;
for j = 1:fDataLen
    force_count = force_count + 1;
    sync_forces(force_count).x_axis = forces(j).x_axis;
    sync_forces(force_count).y_axis = forces(j).y_axis;
    sync_forces(force_count).z_axis = forces(j).z_axis;    
end

ft_count = 1;
force_time = 0;
sync_forces(1).time_diff = force_time;
forceIntervalTime = 0.001;
for k = 2:fDataLen
    ft_count = ft_count + 1;
    force_time = force_time + forceIntervalTime;
    sync_forces(ft_count).time_diff = force_time;
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
    sync_states(hap_count).drill_states = hapticsStates(m).drill_states;
end

dataLen = length(sync_states);
for n = 2:dataLen
    sync_states(n).time_diff = sync_states(n).time_stamp - sync_states(1).time_stamp;
end
sync_states(1).time_diff = 0;

    





