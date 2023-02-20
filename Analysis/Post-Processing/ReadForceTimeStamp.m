function forceTimes = ReadForceTimeStamp(forceTimeStamp)
% ReadForceTimeStamp parses time stamp data recorded from the F/T sensor
% into the struct forceTime
% Load the time stamp file saved by F/T Sensor LabView application into a
% variable using "fopen" function, and then use ReadForceTimeStamp to parse
% data
%
% The struct forceTime saves variables such as hour, minute, second of each
% time stamp line, as well as, the starting and ending time of the data
% recording process. A task duration is also calculated from this data.
%
% e.g. forceTimeStamp = load('filename.txt');
% forceTimes = ReadForceTimeStamp(forceTimeStamp);

timeStamp_count = 0;
while ~feof(forceTimeStamp)
    timeStamp_count = timeStamp_count + 1;
    tline = fgetl(forceTimeStamp);
    tline = strsplit(tline);
    timePeriod = tline(1,3);
    tline = tline(1,2);    
    time = strsplit(tline{1,1},':');
    time = str2double(time);
    forceTimes(timeStamp_count).hour = time(1,1);
    forceTimes(timeStamp_count).minute = time(1,2);
    forceTimes(timeStamp_count).second = time(1,3);
    forceTimes(timeStamp_count).time_period = timePeriod;
end

% Convert to 24 hour format
for i = 1:length(forceTimes)
    if forceTimes(i).time_period == "PM" && forceTimes(i).hour ~= 12
        forceTimes(i).hour = forceTimes(i).hour + 12;
    else 
        forceTimes(i).hour = forceTimes(i).hour;
    end
end

start_hour = forceTimes(1).hour;
start_minute = forceTimes(1).minute;
start_second = forceTimes(1).second;

forceTimes(1).start_hour = start_hour;
forceTimes(1).start_minute = start_minute;
forceTimes(1).start_second = start_second;

end_hour = forceTimes(end).hour;
end_minute = forceTimes(end).minute;
end_second = forceTimes(end).second;

forceTimes(1).end_hour = end_hour;
forceTimes(1).end_minute = end_minute;
forceTimes(1).end_second = end_second;

if forceTimes(1).start_second > forceTimes(1).end_second
    forceTimes(1).end_minute = forceTimes(1).end_minute - 1;
    forceTimes(1).end_second = forceTimes(1).end_second + 60;
end
forceTimes(1).diff_second = forceTimes(1).end_second - forceTimes(1).start_second;

if forceTimes(1).start_minute > forceTimes(1).end_minute
    forceTimes(1).end_hour = forceTimes(1).end_hour - 1;
    forceTimes(1).end_minute = forceTimes(1).end_minute + 60;
end
forceTimes(1).diff_minute = forceTimes(1).end_minute - forceTimes(1).start_minute;
forceTimes(1).diff_hour = forceTimes(1).end_hour - forceTimes(1).start_hour; 

forceTimes(1).duration = (forceTimes(1).diff_hour*3600) + (forceTimes(1).diff_minute*60) + forceTimes(1).diff_second;

fclose(forceTimeStamp);