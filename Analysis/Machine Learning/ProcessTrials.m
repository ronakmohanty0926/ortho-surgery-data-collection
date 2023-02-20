function ProcessTrials(rand_filename, nrand_filename, rand_del_filename, nrand_del_filename)
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

% Setting initial values for the data for easy concatenation
% delta_positions(1,1) = 0;
positions(1,3) = 0;
% delta_speeds(1,1) = 0;
speeds(1,1) = 0;
accelerations(1,1) = 0;
jerks(1,1) = 0;
speed_ratios(1,1) = 0;
% delta_forces(1,1) = 0;
tmp_forces(1,1) = 0;
drillStates(1,1) = 0;
tmp_newtimeForce(1,1) = 0;
timeHaptics(1,1) = 0;

% trainingDataset_number
% For loop to gather training data for each trial
dataLen = length(posFileList);
for i = 1:dataLen
    hapData = posFileList(i);
    hapTime = posTimeFileList(i);
    forceData = forceFileList(i);
    forceTime = forceTimeFileList(i);

    [drillPosition, sync_linearSpeeds, sync_acceleration, sync_jerk, sync_speedRatio, netForces, timeHap, newtimeForce, drill_states] = ProcessTrial(hapData, hapTime, forceData, forceTime);
%     [drillPosition, sync_linearSpeeds, sync_acceleration, sync_jerk, sync_speedRatio, netForces, timeHap, newtimeForce] = ComputeAllParam(device, sync_states, sync_forces);
    %concatenate vertically data for all trials into a single column variable
%     delta_positions = vertcat(delta_positions, delta_position);
%     delta_speeds = vertcat(delta_speeds, delta_speed);
%     delta_forces = vertcat(delta_forces, delta_force);
    drillStates = vertcat(drillStates, drill_states);
    positions = vertcat(positions, drillPosition);
    speeds = vertcat(speeds, sync_linearSpeeds);
    accelerations = vertcat(accelerations, sync_acceleration);
    jerks = vertcat(jerks, sync_jerk);
    speed_ratios = vertcat(speed_ratios, sync_speedRatio);
    tmp_forces = vertcat(tmp_forces, netForces);
    tmp_newtimeForce = vertcat(tmp_newtimeForce, newtimeForce);
    timeHaptics = vertcat(timeHaptics, timeHap);
   
    fclose('all');
%     clear delta_position;
%     clear delta_speed;
%     clear delta_force;
    clear drillPosition; 
    clear sync_linearSpeeds;
    clear sync_acceleration;
    clear sync_jerk;
    clear sync_speedRatio;
    clear netForces;
    clear timeHap;
    clear newtimeForce;
end

forces = tmp_forces(2:(length(positions)+1),1);
timeForces = tmp_newtimeForce(2:(length(positions)+1),1);


%Transfer all columns to a table
AllData = zeros(length(positions),8);
AllData(:,1:3) = positions;
AllData(:,4) = speeds;
AllData(:,5) = accelerations;
AllData(:,6) = jerks;
AllData(:,7) = speed_ratios;
AllData(:,8) = forces;
% AllData(:,8) = timeHaptics;
% AllData(:,10) = timeForces;
AllData(:,9) = drillStates;
AllData(1,:) = [];
writematrix(AllData,nrand_filename);

% [a,b] = size(AllData);
% del_AllData = zeros(a-1, b);
% for k = 2:length(AllData)
%     del_AllData(k-1,1:8) = AllData(k,1:8) - AllData(k-1,1:8);
% end
% del_AllData(:,9) = AllData(1:length(del_AllData),9);
% writematrix(del_AllData,nrand_del_filename);
% 
% %Randomize Data
% [m, n] = size(AllData);
% rand_AllData = zeros(m,n);
% idx = randperm (m); 
% idx = idx';
% for j = 1:length(idx)
%     tmp = idx(j,1);
%     rand_AllData(j,:) = AllData(tmp,:);
% end
% %Write table to file format of your choice
% writematrix(rand_AllData,rand_filename);
% 
% [p, q] = size(del_AllData);
% rand_del_AllData = zeros(p,q);
% index = randperm (p); 
% index = index';
% for s = 1:length(index)
%     tmp = index(s,1);
%     rand_del_AllData(s,:) = del_AllData(tmp,:);
% end
% %Write table to file format of your choice
% writematrix(rand_del_AllData,rand_del_filename);

