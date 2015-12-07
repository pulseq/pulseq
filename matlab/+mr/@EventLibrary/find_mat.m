function [id, found] = find_mat(keys,data,lengths,newData)
    %find_mat Lookup a data structure in the given library.
    %   idx=find_mat(keys,data,lengths,newDat) Return the index of the
    %   newDat in the library. If the data doesn't exist in the library
    %   then the index for the next new entry is returned.
    %
    %   This function is compatible with the MATLAB Coder to generate a mex
    %   file for faster execution.
    %
    %   See also  EventLibrary

    found=false;
    id=0;

    for i=1:length(data)
        if (lengths(i)==length(newData) && norm(data(i).array-newData)<1e-6)
            id=keys(i);
            found=true;
            break;
        end
    end

    if isempty(keys)
        id=1;
    elseif ~found
        id=max(keys)+1;
    end

end
