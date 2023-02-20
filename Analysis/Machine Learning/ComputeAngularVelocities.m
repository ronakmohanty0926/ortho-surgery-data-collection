function angular_velocities = ComputeAngularVelocities(sync_states)
% ComputeAngularVelocities computes angular velocity i.e. rate of change
% of individual gimbal angles per time increment. Read the device state
% struct into this function and an angular velocity struct 
%
% angular_velocities = ComputeAngularVelocities(states)

for i = 2:length(sync_states)
    timeDiff = sync_states(i).time_diff - sync_states(i-1).time_diff;
    gammaDiff = sync_states(i).stylus_gamma - sync_states(i-1).stylus_gamma;
    alphaDiff = sync_states(i).stylus_alpha - sync_states(i-1).stylus_alpha;
    betaDiff = sync_states(i).stylus_beta - sync_states(i-1).stylus_beta;
    angular_velocities(i).alpha = ComputeAngularVelocity(alphaDiff, timeDiff);
    angular_velocities(i).beta = ComputeAngularVelocity(betaDiff, timeDiff);
    angular_velocities(i).gamma = ComputeAngularVelocity(gammaDiff, timeDiff);
end
angular_velocities(1).alpha = 0;
angular_velocities(1).beta = 0;
angular_velocities(1).gamma = 0;

for i = 1:length(angular_velocities)
    temp_alpha(i,1) = angular_velocities(i).alpha;
    temp_beta(i,1) = angular_velocities(i).beta;
    temp_gamma(i,1) = angular_velocities(i).gamma;
end

% temp_alpha = movmean(temp_alpha,125);
% temp_beta = movmean(temp_beta,125);
% temp_gamma = movmean(temp_gamma,125);
temp_alpha = sgolayfilt(temp_alpha,2,125);
temp_beta = sgolayfilt(temp_beta,2,125);
temp_gamma = sgolayfilt(temp_gamma,2,125);

for i = 1:length(angular_velocities)
    angular_velocities(i).alpha = temp_alpha(i,1);
    angular_velocities(i).beta = temp_beta(i,1);
    angular_velocities(i).gamma = temp_gamma(i,1);
end
