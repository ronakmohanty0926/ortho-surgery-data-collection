function [temperature, thermal_times] = ReadTemperatureData(temperatureDataTable)

data_count = 0;
while ~feof(temperatureDataTable)
    data_count = data_count + 1;
    tline = fgetl(temperatureDataTable);
    tline = strsplit(tline);
    timePeriod = tline(1,2);
    temp_time = strsplit(tline{1,1},':');
    temp_time = str2double(temp_time);
    temp_cort1 = tline(1,3);
    temp_cort2 = tline(1,4);
    temp_cort1 = str2double(temp_cort1);
    temp_cort2 = str2double(temp_cort2);
    thermal_times(data_count).hour = temp_time(1,1);
    thermal_times(data_count).minute = temp_time(1,2);
    thermal_times(data_count).second = temp_time(1,3);
    thermal_times(data_count).time_period = timePeriod;
    temperature(data_count).cortex1= temp_cort1;
    temperature(data_count).cortex2= temp_cort2;    
end

% Convert to 24 hour format
for i = 1:length(thermal_times)
    if thermal_times(i).time_period == "PM" && thermal_times(i).hour ~= 12
        thermal_times(i).hour = thermal_times(i).hour + 12;
    else 
        thermal_times(i).hour = thermal_times(i).hour;
    end
end

start_hour = thermal_times(1).hour;
start_minute = thermal_times(1).minute;
start_second = thermal_times(1).second;

thermal_times(1).start_hour = start_hour;
thermal_times(1).start_minute = start_minute;
thermal_times(1).start_second = start_second;

end_hour = thermal_times(end).hour;
end_minute = thermal_times(end).minute;
end_second = thermal_times(end).second;

thermal_times(1).end_hour = end_hour;
thermal_times(1).end_minute = end_minute;
thermal_times(1).end_second = end_second;

if thermal_times(1).start_second > thermal_times(1).end_second
    thermal_times(1).end_minute = thermal_times(1).end_minute - 1;
    thermal_times(1).end_second = thermal_times(1).end_second + 60;
end
thermal_times(1).diff_second = thermal_times(1).end_second - thermal_times(1).start_second;

if thermal_times(1).start_minute > thermal_times(1).end_minute
    thermal_times(1).end_hour = thermal_times(1).end_hour - 1;
    thermal_times(1).end_minute = thermal_times(1).end_minute + 60;
end
thermal_times(1).diff_minute = thermal_times(1).end_minute - thermal_times(1).start_minute;
thermal_times(1).diff_hour = thermal_times(1).end_hour - thermal_times(1).start_hour; 

thermal_times(1).duration = (thermal_times(1).diff_hour*3600) + (thermal_times(1).diff_minute*60) + thermal_times(1).diff_second;


fclose(temperatureDataTable);
