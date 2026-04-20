function out = containers_map_lookup_ik(map,key,fallback)
%CONTAINERS_MAP_LOOKUP this function simulates the missing
%containers.Map.lookup() function using the isKey function
%(supposingly much faster on Octave than exception-try-catch)
    if map.isKey(key)
        out=map(key);
    else
        out=fallback;
    end
end

