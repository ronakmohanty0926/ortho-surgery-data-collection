function [forces, torques] = ReadForceTorqueData(forceDataTable)
% ReadForceTorqueData parses raw Force & Torque data from the F/T sensor \
% into corresponding structs
% Load a file and save it into a variable, the use ReadForceTorque Data to 
% parse data into defined structs: forces and torques
%
% 
% eg. forceDataTable = load('filename.txt')
% [forces, torques] = ReadForceTorqueData(forceDataTable)

dataLen = length(forceDataTable(:,1));

for i = 1:dataLen
    forces(i).x_axis = forceDataTable(i,1);
    forces(i).y_axis = forceDataTable(i,2);
    forces(i).z_axis = forceDataTable(i,3);
    torques(i).x_axis = forceDataTable(i,4);
    torques(i).y_axis = forceDataTable(i,5);
    torques(i).z_axis = forceDataTable(i,6);
end

