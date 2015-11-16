function c=block2events(b)
%BLOCK2EVENTS Convert a block structure to cell array.
%   c=BLOCK2EVENTS(b) Convert the block structure to cell of
%   sequence events.
%
%   If b is already a cell array of events this array is
%   returned unmodified.

c=b;    % Assume b is already a cell array of events
if iscell(b)
    first=b{1}; % Use first element to test for block structure
end
if isfield(first,'rf')
    % Argument is a block structure, copy events to cell array
    % varargin for further processing.
    assert(length(b)==1,'Only a single block structre can be added');
    c={};
    fields=fieldnames(first)';
    for f=fields
        if ~isempty(first.(char(f)))
            c{end+1}=first.(char(f));
        end
    end
elseif iscell(first)
    c=first;
end

end