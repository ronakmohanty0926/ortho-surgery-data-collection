function speedRatio = ComputeSpeedRatio(linearSpeeds)
for i = 2:length(linearSpeeds)
    speedRatio(i,1) = abs(linearSpeeds(i,1) - linearSpeeds(i-1,1));
end
speedRatio(1,1) = 1; 