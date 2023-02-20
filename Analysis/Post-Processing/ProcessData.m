function ProcessData(filename)

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

dataLen = length(posFileList);
% while(1)
% for i = 1:1:dataLen
% %     AllData = [];
% %     count = count + 1;
% %     if(count > 4)
% %         break;
% %     end
%     hapData = posFileList(i);
%     disp(hapData);
%     hapTime = posTimeFileList(i);
%     forceData = forceFileList(i);
%     forceTime = forceTimeFileList(i);
%     [arcLength, signature] =  AnalyzeData(hapData, hapTime, forceData, forceTime);
%     signatures = horzcat(signatures, signature);
%     arcLengths = horzcat(arcLengths, arcLength);
% end
% writematrix(arcLengths, filename,'Sheet',1);  
% writematrix(signatures, filename,'Sheet',2); 
% AllData(:,1) = arcLength;
% AllData(:,2) = signature;

count = 0;
while(1)
    AllData = [];
    count = count + 1;
    if (count > dataLen)
        break;
    end    
    hapData = posFileList(count);
    disp(hapData);
    hapTime = posTimeFileList(count);
    forceData = forceFileList(count);
    forceTime = forceTimeFileList(count);
    [drillPosition, sync_linearSpeeds, sync_netForces, timeHap, timeForce, drill_states] =  AnalyzeData(hapData, hapTime, forceData, forceTime);  
    if (length(timeForce) > length(drillPosition))
        nTimeForce = timeForce(1:length(drillPosition(:,1)),1);
        n_netForce = sync_netForces(1:length(drillPosition(:,1)),1);
    else
        nTimeForce = timeForce;
        n_netForce = sync_netForces;
    end

%         AllData(:,1) = drill_states;
        AllData(:,1:3) = drillPosition;
        AllData(:,4) = timeHap;
%         AllData(:,5) = sync_linearSpeeds;
%         AllData(:,6) = n_netForce;
%         AllData(:,7) = timeHap;
%         AllData(:,8) = nTimeForce;
        
%     nCount = 0;    
%     for i = 1:length(AllData)
%         if (AllData(i,1) == 1 || AllData(i,1) == 2 || AllData(i,1) == 3)
%             nCount = nCount + 1;
%             newData(nCount,:) = AllData(i,:);
%         end
%     end
            
         
    writematrix(AllData, filename,'Sheet',count);    
%     close('all');
end



