function netForces = ComputeNetForces(sync_forces)
% ComputeNetForces computes the net forces for a drilling trial. The forces
% struct is given as an input with forces along individual axes and the net
% forces are computed using ComputeNetForces.
%
% e.g. netForces = ComputeNetForces(forces)

dataLen = length(sync_forces);

for i = 1:dataLen
    netForces(i,1) = ComputeNetForce(sync_forces(i));
end
% netForces = movmean(netForces,125);
% netForces = sgolayfilt(netForces,2,125);