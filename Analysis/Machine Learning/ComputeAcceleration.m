function acceleration = ComputeAcceleration(linearSpeeds, timeHap)

for i = 2:length(linearSpeeds)
    delta_speed = linearSpeeds(i,1) - linearSpeeds(i-1,1);
    delta_time = timeHap(i,1) - timeHap(i-1,1);
    acceleration(i,1) = abs(delta_speed./delta_time);
end
acceleration(1,1) = 0; 
    
    