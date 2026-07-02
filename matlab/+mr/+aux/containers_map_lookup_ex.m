function out = containers_map_lookup_ex(map,key,fallback)
%CONTAINERS_MAP_LOOKUP this function simulates the missing
%containers.Map.lookup() function using try-catch clause
    try
        out=map(key);
    catch
        out=fallback;
    end
end

