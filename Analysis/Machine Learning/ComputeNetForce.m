function netForce = ComputeNetForce(force)
% ComputeNetForce computes the net force at given time frame. The force
% struct is given as an input with forces along individual axes and the net
% force is computed using ComputeNetForce.
%
% e.g. netForce = ComputeNetForce(force)

sq_netForce = (force.x_axis^2) + (force.y_axis^2) + (force.z_axis^2);

netForce = sqrt(sq_netForce);
