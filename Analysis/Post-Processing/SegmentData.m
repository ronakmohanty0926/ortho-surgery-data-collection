function labelled_data = SegmentData(drill_states, sync_netForces, sync_linearSpeeds, drillt_coordinates, timeHap, timeForce)

count = 0;
skipFrame = 3;
idx = 1;
netForces(1,1) = sync_netForces(1,1);
newtimeForce(1,1) = timeForce(1,1);

while count < (length(sync_netForces) - skipFrame)
    count = count + skipFrame;
    idx = idx + 1;
    netForces(idx,1) = sync_netForces(count,1);
    newtimeForce(idx,1) = timeForce(count,1);
end

idx_away = 0;
idx_cort1 = 0;
idx_cort2 = 0;
idx_retract = 0;

for i = 1:length(drill_states)
    if drill_states(i,1) == 0
        idx_away = idx_away + 1;
        position_away(idx_away,1) = drillt_coordinates(i,1);
        timeHap_away(idx_away,1) = timeHap(i,1);
        speed_away(idx_away,1) = sync_linearSpeeds(i,1);       
        force_away(idx_away,1) = netForces(i,1);
        timeForce_away(idx_away,1) = newtimeForce(i,1);
    elseif drill_states(i,1) == 1
        idx_cort1 = idx_cort1 + 1;
        labelled_data(idx_cort1).position_cort1 = drillt_coordinates(i,1);
        labelled_data(idx_cort1).timeHap_cort1 = timeHap(i,1);
        labelled_data(idx_cort1).speed_cort1 = sync_linearSpeeds(i,1);
        labelled_data(idx_cort1).force_cort1 = netForces(i,1);
        labelled_data(idx_cort1).timeForce_cort1 = newtimeForce(i,1);
    elseif drill_states(i,1) == 2
        idx_cort2 = idx_cort2 + 1;
        labelled_data(idx_cort2).position_cort2 = drillt_coordinates(i,1);
        labelled_data(idx_cort2).timeHap_cort2 = timeHap(i,1);
        labelled_data(idx_cort2).speed_cort2 = sync_linearSpeeds(i,1);
        labelled_data(idx_cort2).force_cort2 = netForces(i,1);
        labelled_data(idx_cort2).timeForce_cort2 = newtimeForce(i,1);
    elseif drill_states(i,1) == 3
        idx_retract = idx_retract + 1;
        labelled_data(idx_retract).position_retract = drillt_coordinates(i,1);
        labelled_data(idx_retract).timeHap_retract = timeHap(i,1);
        labelled_data(idx_retract).speed_retract = sync_linearSpeeds(i,1);
        labelled_data(idx_retract).force_retract = netForces(i,1);
        labelled_data(idx_retract).timeForce_retract = newtimeForce(i,1);
    end
end

idx_approach = 0;
idx_rest = 0;
for j = 1:length(position_away)
    if timeHap_away(j,1) < labelled_data(1).timeHap_cort1
        idx_approach = idx_approach + 1;
        labelled_data(idx_approach).position_approach = position_away(j,1);
        labelled_data(idx_approach).timeHap_approach = timeHap_away(j,1);
        labelled_data(idx_approach).speed_approach = speed_away(j,1);
    elseif timeHap_away(j,1) > labelled_data(1).timeHap_cort1
        idx_rest = idx_rest + 1;
        labelled_data(idx_rest).position_rest = position_away(j,1);
        labelled_data(idx_rest).timeHap_rest = timeHap_away(j,1);
        labelled_data(idx_rest).speed_rest = speed_away(j,1);
    end
end

idx_fapproach = 0;
idx_frest = 0;
for k = 1:length(force_away)
    if timeForce_away(k,1) < labelled_data(1).timeForce_cort1
        idx_fapproach = idx_fapproach + 1;
        labelled_data(idx_fapproach).force_approach = force_away(k,1);
        labelled_data(idx_fapproach).timeForce_approach = timeForce_away(k,1);
    elseif timeForce_away(k,1) > labelled_data(1).timeForce_cort1
        idx_frest = idx_frest + 1;
        labelled_data(idx_frest).force_rest = force_away(k,1);
        labelled_data(idx_frest).timeForce_rest = timeForce_away(k,1);
    end
end



