function angular_velocity = ComputeAngularVelocity(angleDiff, timeDiff)
% ComputeAngularVelocity computes the angular velocity for the rate of
% change of angle per time increment from its previous time frame. The
% inputs provided are the change in angle and time increment resulting
% which is further processed by ComputeAngularVelocity to compute the
% angular velocity.
%
%e.g. angular_velocity = ComputeAngularVelocity(angleDiff, timeDiff)

angular_velocity = angleDiff./timeDiff;



