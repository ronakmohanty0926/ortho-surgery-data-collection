function linearSpeeds = ComputeLinearSpeeds(position, timeHap)

dataLen = length(position);



% for i = 301:dataLen
for j = 2:dataLen     
    movX = position(j,1) - position(j-1,1);
    movY = position(j,2) - position(j-1,2);
    movZ = position(j,3) - position(j-1,3);
    
%     movX = posX(i,1) - posX(i-300,1);
%     movY = posY(i,1) - posY(i-300,1);
%     movZ = posZ(i,1) - posZ(i-300,1);
    
%     timeDiff = sync_states(i).time_diff - sync_states(i-300).time_diff;
    timeDiff = timeHap(j,1) - timeHap(j-1,1);
    
    velX = movX./timeDiff;
    velY = movY./timeDiff;
    velZ = movZ./timeDiff;
    
    sq_linearSpeed = (velX^2) + (velY^2) + (velZ^2);

    linearSpeeds(j,1) = abs(sqrt(sq_linearSpeed));  
end
linearSpeeds(1,1) = 0;
% linearSpeeds(251,1) = 0;
% linearSpeeds = movmean(linearSpeeds,125);
% linearSpeeds = sgolayfilt(linearSpeeds,2,125);


    





