function DeviceUI(device, state, option)
% '[h, coordinates] = DisplayDeviceState(device, state, option)' displays 
% the Geomagic device for a single state at a given instance of time.
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
% 'option' can be 'basic' or 'fancy' and determines the visualization type.
% The 'basic' option renders a simple line-schematic of the device. The
% 'fancy' option renders a 3D version (spheres and cylinders) of the device.
%
% Example:
%
% device.link_length1 = 135;
% device.link_length2 = 135;
% device.tip_length = 35;
% device.top_length = 55;
% device.base_height = 160;
% 
% state.base_angle = pi/6;
% state.angle_link01 = pi/4;
% state.angle_link12 = -pi/6;
% state.stylus_alpha = 0;
% state.stylus_beta = 0;
% state.stylus_gamma = 0;
% 
% [h, coordinates] = DisplayDeviceState(device, state, 'fancy');

h = figure('Units','normalized','Position',[0 0 1 1]);
view(3)
camup([0 1 0])
hold on;
Rng = device.link_length1+device.link_length2;
axis([-1.2*Rng 1.2*Rng -1.2*device.base_height 1.2*Rng -100 1.2*Rng]);
set(gca, 'visible', 'off');
set(gcf,'color','w');
campos([1*Rng 1*Rng-0.5*device.base_height 1*Rng])
daspect([1 1 1]);
light('Position',[0.5*Rng Rng-0.5*device.base_height 0.5*Rng],'Style','infinite');

coordinates = ComputeCoordinatesForState(device, state);

% Display axes
line([0 0.5*Rng], [0 0], [0 0], 'Color','red', 'LineWidth', 4);
line([0 0], [0 Rng], [0 0], 'Color','green', 'LineWidth', 4);
line([0 0], [0 0], [0 Rng], 'Color','blue', 'LineWidth', 4);

% Display floor
[floor_x, floor_z] = meshgrid(-1.2*Rng:0.6*Rng:1.2*Rng); % Generate x and z data
floor_y = -device.base_height*ones(size(floor_x, 1)); % Generate y data
surf(floor_x, floor_y, floor_z,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7)

if(option == 'basic')
    line([0 coordinates.elbow(1)], [0 coordinates.elbow(2)], [0 coordinates.elbow(3)], 'Color','black', 'LineWidth', 4);
    line([coordinates.elbow(1) coordinates.wrist(1)], [coordinates.elbow(2) coordinates.wrist(2)], [coordinates.elbow(3) coordinates.wrist(3)], 'Color','black', 'LineWidth', 4);

    line([coordinates.wrist(1) coordinates.tip_stylus(1)], [coordinates.wrist(2) coordinates.tip_stylus(2)], [coordinates.wrist(3) coordinates.tip_stylus(3)], 'Color',[1.0 0.55 0.0], 'LineWidth', 4);
    line([coordinates.wrist(1) coordinates.top_stylus(1)], [coordinates.wrist(2) coordinates.top_stylus(2)], [coordinates.wrist(3) coordinates.top_stylus(3)], 'Color',[1.0 0.55 0.0], 'LineWidth', 4);
elseif(option == 'fancy')        
    % Display Joints using Shperes
    [X, Y, Z] = sphere(20);
%    surf(50.*X(:,11:end),50.*Y(:,11:end),50.*Z(:,11:end),'FaceColor',[0.2 0.2 0.6],'EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.85);
    surf(50.*X,50.*Y,50.*Z,'FaceColor',[0.2 0.2 0.6],'EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.85);
    surf(20.*X + coordinates.elbow(1),20.*Y + coordinates.elbow(2),20.*Z + coordinates.elbow(3),'FaceColor',[0.2 0.2 0.6],'EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.85);
    surf(10.*X + coordinates.wrist(1),10.*Y + coordinates.wrist(2),10.*Z + coordinates.wrist(3),'FaceColor',[0.2 0.2 0.6],'EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.85);
    
    % Display Base
    [Xb, Yb, Zb] = cylinder2P([70, 65, 50, 35, 15], 100, 5, [0 -device.base_height 0], [0 0 0]);
    surf(Xb,Yb,Zb,'FaceColor',[0.2 0.2 0.2],'FaceLighting','gouraud','EdgeColor','none','AmbientStrength',1.0);
    
    % Display Links using Cylinders
    sc = Rng/7;
    [X1, Y1, Z1] = cylinder2P(sc.*[0.25,0.25,0.20], 25, 3, [0 0 0], coordinates.elbow');
    [X2, Y2, Z2] = cylinder2P(sc.*[0.20,0.20,0.15], 25, 3, coordinates.elbow', coordinates.wrist');
    [X3, Y3, Z3] = cylinder2P(sc.*[0.15,0.15,0.15], 25, 3, coordinates.wrist', coordinates.top_stylus');
    [X4, Y4, Z4] = cylinder2P(sc.*[0.15,0.15,0.05], 25, 3, coordinates.wrist', coordinates.tip_stylus');
    
    surf(X1,Y1,Z1,'FaceColor',[0.5 0.5 0.5],'FaceLighting','gouraud','EdgeColor','none','AmbientStrength',0.75);
    surf(X2,Y2,Z2,'FaceColor',[0.5 0.5 0.5],'FaceLighting','gouraud','EdgeColor','none','AmbientStrength',0.75);
    surf(X3,Y3,Z3,'FaceColor',[1.0 0.55 0.0],'FaceLighting','gouraud','EdgeColor','none','AmbientStrength',0.85);
    surf(X4,Y4,Z4,'FaceColor',[1.0 0.55 0.0],'FaceLighting','gouraud','EdgeColor','none','AmbientStrength',0.85);
    
    % light('Position',[Rng Rng-0.5*device.base_height Rng],'Style','local');
end