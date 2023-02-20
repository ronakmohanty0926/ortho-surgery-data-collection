function linearSpeeds = ComputeLinearSpeeds(position, timeHap)

dataLen = length(position);



for i = 301:dataLen
% for j = 2:dataLen     
%     movX = position(j,1) - position(j-1,1);
%     movY = position(j,2) - position(j-1,2);
%     movZ = position(j,3) - position(j-1,3);
    
    movX = position(i,1) - position(i-1,1);
    movY = position(i,2) - position(i-1,2);
    movZ = position(i,3) - position(i-1,3);
    
    timeDiff = timeHap(i) - timeHap(i-300);
%     timeDiff = timeHap(j,1) - timeHap(j-1,1);
    
    velX = movX./timeDiff;
    velY = movY./timeDiff;
    velZ = movZ./timeDiff;
    
    sq_linearSpeed = (velX^2) + (velY^2) + (velZ^2);

    linearSpeeds(i,1) = abs(sqrt(sq_linearSpeed));  
%     linearSpeeds(j,1) = abs(sqrt(sq_linearSpeed));  
end
linearSpeeds(1,1) = 0;
% linearSpeeds(251,1) = 0;
% linearSpeeds = movmean(linearSpeeds,125);
% linearSpeeds = sgolayfilt(linearSpeeds,2,125);


    





