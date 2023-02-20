% ProcessTrials('Train_Test_Force.xlsx', 'Train_Test_Force_nRand.xlsx', 'Train_Test_Del_Force.xlsx', 'Train_Test_Del_Force_nRand.xlsx')
% ProcessTrials('Test_Force.xlsx', 'Test_Force_nRand.xlsx', 'Test_Del_Force.xlsx', 'Test_Del_Force_nRand.xlsx')
% AnalyzeData('Log_Kellam_6.txt', 'Time_Kellam_6.txt','Force_6.txt', 'Time_6.txt')
%[delta_position, delta_speed, delta_force, drill_states] = ProcessTrial('Log_Kellam_6.txt', 'Time_Kellam_6.txt','Force_6.txt', 'Time_6.txt');
ComputeDimension('OB_new.mat')
