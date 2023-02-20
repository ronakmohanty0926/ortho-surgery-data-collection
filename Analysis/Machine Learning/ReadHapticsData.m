function [device, hapticsStates] = ReadHapticsData(hapticsDataTable)
% ReadHapticsData parses haptics raw data into corresponding structs
% Load a file and save it into a variable, the use ReadHapticsData to parse data
% into defined structs: device and states
%
% device struct stores device measurement variables
% hapticsStates struct stores the RAW encoder data given by the haptics API
%
% eg. hapticsDataTable = load('filename.txt')
% [device, hapticsStates] = ReadHapticsData(hapticsDataTable)

frameSize = length(hapticsDataTable(:,1));

device.link_length1 = 135;
device.link_length2 = 130;
device.tip_length = 40;
device.top_length = 135;
device.base_height = 163;
device.button1_distance = 15;
device.button1_height = 5;
device.button2_distance = 30;
device.button2_height = 5;

device.drilltip_offseth = 163; % drill tip horizontal offset from the gimbal
device.drilltip_offsetv = 45.62; % drill tip vertical offset from the gimbal
device.drilltop_offseth = 115.0; % drill top horizontal offset from the gimbal
device.drilltop_offsetv = 45.62; % drill top vertical offset from the gimbal
device.holder_offseth = 75; % drill holder horizontal offset from the gimbal
device.holder_offsetv = 45.62; % drill holder horizontal offset from the gimbal
device.drillfixture_offseth = 75.0;
device.drillbase_length = 220;
device.forcesensor_height = 131;

for i = 1:frameSize
    
    hapticsStates(i).base_angle = -1*hapticsDataTable(i,18);
    hapticsStates(i).angle_link01 = -1.*(hapticsDataTable(i,19));
    hapticsStates(i).angle_link12 = (pi./2) - (hapticsDataTable(i,20) - hapticsDataTable(i,19)); %this angle is recorded wrt global coord system
    
    hapticsStates(i).stylus_gamma = -1.*hapticsDataTable(i,21);
    hapticsStates(i).stylus_alpha = -1.*hapticsDataTable(i,22);
    hapticsStates(i).stylus_beta = -1.*hapticsDataTable(i,23);    
    
    hapticsStates(i).time_stamp = hapticsDataTable(i,2);
    
    hapticsStates(i).drill_states = hapticsDataTable(i,1);
end