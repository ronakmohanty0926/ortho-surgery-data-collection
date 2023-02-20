function coordinates = ComputeCoordinatesForState(device, state)
% 'coordinates = ComputecoordinatesFromState(device, state)' computes the position 
% vectors of each of the links of the Geomagic Device for the given state 
% variables.
% 
% 'device' is a structure that contains the parameters (link lengths) of the 
% Geomagic device, namely, 'link_length1', 'link_length2', 'tip_length',
% and 'top_length'. Each of these variables are 1x1 vectors. All lengths
% are in millimeters.
%
% 'state' is a structure that contains the all angles of the Geomagic device
% at a given instance of time. The structure contains six angle parameters, 
% namely, 'base_angle', 'angle_link01', 'angle_link12', 'stylus_alpha', 'stylus_beta',
% 'stylus_gamma', and 'time_stamp'. Each of these variables are 1x1
% vectors. All angles are in radians and time-stamp is in seconds.
%
% 'coordinates' is a structure that contains the coordinates of the joints at
% the end of each link, namely, 'elbow', 'wrist', 'tip_stylus', 
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
% state.base_angle = pi/3;
% state.angle_link01 = pi/4;
% state.angle_link12 = -pi/6;
% state.stylus_alpha = pi/4;
% state.stylus_beta = pi/2;
% state.stylus_gamma = pi/10;
% state.timestamp = -1; % this will assume static data
%
% coordinates = ComputecoordinatesFromState(device, state);

base_rotation = BasicRotationMatrix('Y', state.base_angle);

coordinates.elbow = base_rotation*[0;...
    -device.link_length1*sin(state.angle_link01);...
    device.link_length1*cos(state.angle_link01)];

coordinates.wrist = base_rotation*[0;...
    -device.link_length1*sin(state.angle_link01)-device.link_length2*sin(state.angle_link01+state.angle_link12);...
    device.link_length1*cos(state.angle_link01)+device.link_length2*cos(state.angle_link01+state.angle_link12)];

R_alpha = BasicRotationMatrix('X', state.stylus_alpha); % Pitch Angle
R_beta = BasicRotationMatrix('Y', state.stylus_beta); % Yaw Angle
R_gamma = BasicRotationMatrix('Z', state.stylus_gamma); % Roll Angle

G = R_gamma*R_alpha*R_beta; % Roll-Pitch-Yaw Rotation Matrix

R_12 = BasicRotationMatrix('X', state.angle_link12);
p_12 = [0;0;device.link_length2];
x_12 = R_12*(p_12 + G(:,1));
y_12 = R_12*(p_12 + G(:,2));
z_12 = R_12*(p_12 + G(:,3));

R_01 = BasicRotationMatrix('X', state.angle_link01);
p_01 = [0;0;device.link_length1];
x_02 = base_rotation*R_01*(p_01 + x_12);
y_02 = base_rotation*R_01*(p_01 + y_12);
z_02 = base_rotation*R_01*(p_01 + z_12);

diff_x = x_02-coordinates.wrist;
diff_y = y_02-coordinates.wrist;
diff_z = z_02-coordinates.wrist;

coordinates.end_frame = [diff_x./norm(diff_x), diff_y./norm(diff_y), diff_z./norm(diff_z)];

coordinates.tip_stylus = coordinates.wrist - (device.tip_length*coordinates.end_frame(:,2));
coordinates.top_stylus = coordinates.wrist + (device.top_length*coordinates.end_frame(:,2));

coordinates.button1_base = coordinates.wrist + (device.button1_distance*coordinates.end_frame(:,2));
coordinates.button2_base = coordinates.wrist + (device.button2_distance*coordinates.end_frame(:,2));
coordinates.button1 = coordinates.button1_base - (device.button1_height*coordinates.end_frame(:,3));
coordinates.button2 = coordinates.button2_base - (device.button2_height*coordinates.end_frame(:,3));

% coordinates.drill_offseth = coordinates.wrist + (device.drilloffseth*coordinates.end_frame(:,2));
% coordinates.drill_offsetv = coordinates.drill_offseth + (device.drilloffsetv*coordinates.end_frame(:,3));
coordinates.drilltop = coordinates.wrist + (device.drilltop_offseth*coordinates.end_frame(:,2)) + (device.drilltop_offsetv*coordinates.end_frame(:,3));
coordinates.drilltip = coordinates.wrist - (device.drilltip_offseth*coordinates.end_frame(:,2)) + (device.drilltip_offsetv*coordinates.end_frame(:,3));
coordinates.drillfix = coordinates.wrist + (device.drilltip_offsetv*coordinates.end_frame(:,3)) - (device.drillfixture_offseth*coordinates.end_frame(:,2));
coordinates.drillholdertop = coordinates.wrist + (device.holder_offseth*coordinates.end_frame(:,2));
coordinates.drillholderbase = coordinates.drillholdertop + (device.holder_offsetv*coordinates.end_frame(:,3));
coordinates.drillbase = coordinates.drilltop + (device.drillbase_length*coordinates.end_frame(:,3));