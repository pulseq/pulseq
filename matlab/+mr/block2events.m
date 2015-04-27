function c=block2events(b)
%BLOCK2EVENTS Convert a block structure to cell array.
%   c=BLOCK2EVENTS(b) Convert the block structure to cell of
%   sequence events.
%
%   If b is already a cell array of events this array is
%   returned unmodified.

c=b;    % Assume b is already a cell array
if iscell(b)
    b=b{1}; % Check if first element is block structure
end
if isfield(b,'rf')
    % Argument is a block structure, copy events to cell array
    % varargin for further processing.
    s=b;
    c={};
    fields=fieldnames(s)';
    for f=fields
        if ~isempty(s.(char(f)))
            c{end+1}=s.(char(f));
        end
    end
end

end