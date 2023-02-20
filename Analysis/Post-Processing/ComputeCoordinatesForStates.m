function coordinates_series = ComputeCoordinatesForStates(device, states)
% 'coordinates_series = ComputeCoordinatesForStates(device, states)' generates
% the coordinates of the Geomagic device for a time-series of device states.
% 
% 'device' is a structure that contains the parameters (link lengths) of  
% the Geomagic device, namely, 'link_length1', 'link_length2', 'tip_length',
% and 'top_length'. Each of these variables are 1x1 vectors. All lengths
% are in millimeters.
%
% 'states' is a structure-array that contains the all angles of the 
% Geomagic device as a timeseries. Each element 'states(i)' contains six  
% angle parameters, namely, 'base_angle', 'angle_link01', 'angle_link12', 
% 'stylus_alpha', 'stylus_beta', 'stylus_gamma', and 'time_stamp'. Each of 
% these variables are 1x1 vectors. All angles are in radians and 
% time-stamps in seconds.
% 
% 'coordinates_series' is a structure-array that contains the time-series of 
% coordinates of the joints at the end of each link. Each element
% 'coordinates_series(i)' contains 'elbow', 'wrist', 'tip_stylus', 
% 'top_stylus', and 'end_frame'. Each of these variables are 3x1 vectors 
% except 'end_frame' which is a 3x3 matrix that represents the coordinate 
% frame of the end-effector. All coordinates are in millimeters.
%
% Example:
%
% device.link_length1 = 135;
% device.link_length2 = 130;
% device.tip_length = 35;
% device.top_length = 55;
% 
% ang_base = [zeros(1,50) pi/100:pi/100:pi/4];
% ang_var1 = [-pi/100:-pi/100:-pi/4 (-pi/4).*ones(1,50)];
% ang_var2 = [zeros(1,25) pi/100:pi/100:pi/4 (pi/4).*ones(1,25)];
% time_stamps = 0:0.001:0.0749;
% 
% for i = 1:75
%     states(i).base_angle = ang_base(i);
%     states(i).angle_link01 = ang_var1(i);
%     states(i).angle_link12 = ang_var2(i);
%     states(i).stylus_alpha = 0;
%     states(i).stylus_beta = 0;
%     states(i).stylus_gamma = 0;
%     states(i).time_stamp = time_stamps(i);
% end
% 
% coordinates_series = ComputeCoordinatesForStates(device, states);

for i = 1:length(states)
    coordinates_series(i) = ComputeCoordinatesForState(device, states(i));
end