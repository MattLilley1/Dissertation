% This script is a basic interpolation search algorithm used to assign the
% correct stimulus value with the correct time value
% takes 3 parameters:
%   timeArr : array of time value intervals, these are uniformly
%       distributed to represent the 24 hours in the day.
%   target : value to find.

% index can be defined here but it is used in fitting_uv_efff.m
function index = interpolationSearch(timeArr, target)
    low = 1;
    high = length(timeArr);
    
    while low <= high && target >= timeArr(low) && target <= timeArr(high)
        if low == high
            if timeArr(low) == target
                index = low;
                return;
            else
                index = -1;
                return;
            end
        end
        
        pos = low + floor((((high - low) / (timeArr(high) - timeArr(low))) * (target - timeArr(low))));
        
        if timeArr(pos) == target
            index = pos;
            return;
        elseif timeArr(pos) < target
            low = pos + 1;
        else
            high = pos - 1;
        end
    end
    
    index = -1; % Target not found
end