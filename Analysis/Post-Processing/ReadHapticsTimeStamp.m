function hapticsTimes = ReadHapticsTimeStamp(hapticsTimeStamp)
% ReadHapticsTimeStamp parses time stamp data recorded from the position
% sensing device into the struct hapticsTimes
% Load the time stamp file saved by OpenHaptics API into a variable using 
% "fopen" function. Later input this variable and use ReadHapticsTimeStamp 
% to parse data
%
% The struct hapticsTimes saves variables such as the starting and ending time of the data
% recording process. A task duration is also calculated from this data.
%
% e.g. hapticsTimeStamp = load('filename.txt');
% hapticsTimes = ReadHapticsTimeStamp(hapticsTimeStamp);

timeStamp_count = 0;
while ~feof(hapticsTimeStamp)
    timeStamp_count = timeStamp_count + 1;
    tline = fgetl(hapticsTimeStamp);
    if (tline == "")
    {};  
    else
        tline = strsplit(tline);
        tline = tline(1,4);
        time = strsplit(tline{1,1},':');
        time = str2double(time);
        hapticsTimes(timeStamp_count).hour = time(1,1);
        hapticsTimes(timeStamp_count).minute = time(1,2);
        hapticsTimes(timeStamp_count).second = time(1,3);  
    end
end

start_hour = hapticsTimes(1).hour;
start_minute = hapticsTimes(1).minute;
start_second = hapticsTimes(1).second;

hapticsTimes(1).start_hour = start_hour;
hapticsTimes(1).start_minute = start_minute;
hapticsTimes(1).start_second = start_second;

end_hour = hapticsTimes(end).hour;
end_minute = hapticsTimes(end).minute;
end_second = hapticsTimes(end).second;

hapticsTimes(1).end_hour = end_hour;
hapticsTimes(1).end_minute = end_minute;
hapticsTimes(1).end_second = end_second;

if hapticsTimes(1).start_second > hapticsTimes(1).end_second
    hapticsTimes(1).end_minute = hapticsTimes(1).end_minute - 1;
    hapticsTimes(1).end_second = hapticsTimes(1).end_second + 60;
end
hapticsTimes(1).diff_second = hapticsTimes(1).end_second - hapticsTimes(1).start_second;

if hapticsTimes(1).start_minute > hapticsTimes(1).end_minute
    hapticsTimes(1).end_hour = hapticsTimes(1).end_hour - 1;
    hapticsTimes(1).end_minute = hapticsTimes(1).end_minute + 60;
end
hapticsTimes(1).diff_minute = hapticsTimes(1).end_minute - hapticsTimes(1).start_minute;
hapticsTimes(1).diff_hour = hapticsTimes(1).end_hour - hapticsTimes(1).start_hour; 

hapticsTimes(1).duration = (hapticsTimes(1).diff_hour*3600) + (hapticsTimes(1).diff_minute*60) + hapticsTimes(1).diff_second;

fclose(hapticsTimeStamp);
