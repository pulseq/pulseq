function c = block2events(b)
%BLOCK2EVENTS Convert a block structure to cell array.
%   c=BLOCK2EVENTS(b) Convert the block structure to cell of
%   sequence events.
%
%   If b is already a cell array of events this array is
%   returned unmodified.

% strip away 1x1 cell wrapper(s) -- otherwise it conflicts with adding ready-made blocks
while iscell(b) && 1==length(b) && iscell(b{1})
    b=b{1};
end

c = b;    % Assume b is already a cell array of events
if iscell(b)
    first = b{1}; % Use first element to test for block structure
else 
    first = b;
end
if isfield(first, 'rf')
    % Argument is a block structure, copy events to cell array
    % varargin for further processing.
    assert(length(b) == 1, 'Only a single block structure can be added');
%    fields = fieldnames(first)';
%%%%%%%
%     c = {};
%     for f = fields
%         if ~isempty(first.(char(f)))
%             c{end+1} = first.(char(f));
%         end
%     end
%%
%     c = cell(1,length(fields));
%     for i = 1:length(fields)
%         c{i} = first.(char(fields(i)));
%     end
    c=struct2cell(first);
    %c(cellfun(@isempty,c))=[];
    c(cellfun('isempty',c))=[];
    %c=c(~cellfun('isempty',c));
    
    %c(cellfun('isnumeric',c))=[]; % remove blockDuration 

elseif iscell(first)
    c = first;
end

end