function jerk = ComputeJerk(acceleration, timeHap)

for i = 2:length(acceleration)
    delta_acc = acceleration(i,1) - acceleration(i-1,1);
    delta_time = timeHap(i,1) - timeHap(i-1,1);
    jerk(i,1) = abs(delta_acc./delta_time);
end
jerk(1,1) = 0; 