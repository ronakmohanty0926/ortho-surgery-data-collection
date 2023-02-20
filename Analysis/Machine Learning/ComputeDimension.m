function ComputeDimension(filename)

% function ProcessTrials(rand_filename, nrand_filename, rand_del_filename, nrand_del_filename)
%'ProcessTrials(trainingDataset_number, filename)' parses all trial data
% and outputs a matrix as a CSV, EXCEL, TEXT as per choice
%
% used for preparing the training data out of all data files
%
% example: ProcessTrials(5, 'test.xlsx') 

% Looks for the identifier for Haptic Position files and stores the name as
% a string list
tmp_posFileList = dir(fullfile(pwd,'*PO_*'));
tmp_posFileList1 = {tmp_posFileList.name};
posFileList = string(tmp_posFileList1);
posFileList = posFileList';

% Looks for the identifier for Haptic Time files and stores the name as
% a string list
tmp_posTimeFileList = dir(fullfile(pwd,'*PT_*'));
tmp_posTimeFileList = {tmp_posTimeFileList.name};
posTimeFileList = string(tmp_posTimeFileList);
posTimeFileList = posTimeFileList';

% Looks for the identifier for Force files and stores the name as
% a string list
tmp_forceFileList = dir(fullfile(pwd,'*FO_*'));
tmp_forceFileList = {tmp_forceFileList.name};
forceFileList = string(tmp_forceFileList);
forceFileList = forceFileList';

% Looks for the identifier for Force Time files and stores the name as
% a string list
tmp_forceTimeFileList = dir(fullfile(pwd,'*FT_*'));
tmp_forceTimeFileList = {tmp_forceTimeFileList.name};
forceTimeFileList = string(tmp_forceTimeFileList);
forceTimeFileList = forceTimeFileList';

dim1 = zeros(1,1);
dim2 = zeros(1,1);
dataLen = length(posFileList);
count = 0;
while(1)
    count = count + 1;
    if (count > dataLen)
        break;
    end
    hapData = posFileList(count);
    hapTime = posTimeFileList(count);
    forceData = forceFileList(count);
    forceTime = forceTimeFileList(count);


    [drillPosition, sync_linearSpeeds, sync_acceleration, sync_jerk, sync_speedRatio, netForces, timeHap, newtimeForce, drill_states] = ProcessTrial(hapData, hapTime, forceData, forceTime);
%     [drillPosition, sync_linearSpeeds, sync_acceleration, sync_jerk, sync_speedRatio, netForces, timeHap, newtimeForce] = ComputeAllParam(device, sync_states, sync_forces);
    %concatenate vertically data for all trials into a single column variable
%     delta_positions = vertcat(delta_positions, delta_position);
%     delta_speeds = vertcat(delta_speeds, delta_speed);
%     delta_forces = vertcat(delta_forces, delta_force);
%     drillStates = vertcat(drillStates, drill_states);
%     positions = vertcat(positions, drillPosition);
%     speeds = vertcat(speeds, sync_linearSpeeds);
%     accelerations = vertcat(accelerations, sync_acceleration);
%     jerks = vertcat(jerks, sync_jerk);
%     speed_ratios = vertcat(speed_ratios, sync_speedRatio);
%     tmp_forces = vertcat(tmp_forces, netForces);
%     tmp_newtimeForce = vertcat(tmp_newtimeForce, newtimeForce);
%     timeHaptics = vertcat(timeHaptics, timeHap);

AllData = zeros(length(drillPosition),1);
AllData(:,1) = drill_states;
AllData(:,2:4) = drillPosition;
cortex1_count = 0;
cortex2_count = 1;

for i = 1:length(drillPosition)
    if (drill_states(i) == 1)
        cortex1_count = cortex1_count + 1;
        cortex1(cortex1_count,1) = drillPosition(i,2);
    elseif (drill_states (i) == 2)
        cortex2_count = cortex2_count + 1;
        cortex2(cortex2_count,1) = drillPosition(i,2);   
    end
end

cort1_dim = abs(cortex1(end,1) - cortex1(1,1));
cort2_dim = abs(cortex2(end,1) - cortex2(1,1));

dim1 = vertcat(dim1,cort1_dim);
dim2 = vertcat(dim2,cort2_dim);
% writematrix(AllData,filename,'Sheet',count);   
fclose('all');
%     clear delta_position;
%     clear delta_speed;
%     clear delta_force;
    clear drill_states;
    clear drillPosition; 
    clear sync_linearSpeeds;
%     clear sync_acceleration;
%     clear sync_jerk;
%     clear sync_speedRatio;
    clear netForces;
    clear forces;
    clear timeHap;
    clear newtimeForce;
end

save(filename, 'dim1', 'dim2');