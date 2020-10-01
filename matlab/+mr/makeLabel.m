function out = makeLabel(type, label, value)
%makeLabel Create a ADC Label.
%   label=makeLabel() Create a Label event for ADC line
%         Possible values for type are: 'SET','INC'.
%         Label may be a counter of a flag and should be one of
%         (counters) 'SLC','SEG','REP','AVG','SET','ECO','PHS','LIN','PAR', 
%         (flags)    'NAV','REV','SMS'.
%         Value: numeric value of the parameter (increment may be negative)  
%                or true/false for a flag
%
%   See also  Sequence.addBlock, mr.getSupportedLabels


supported_labels=mr.getSupportedLabels();
if nargin~=3
    error('makeLabel:invalidArguments','Must supply exactly 3 parameters');
end
if ~any(ismember(supported_labels,label))
    error('makeLabel:invalidArguments','Must supply a valid label');
end
if ~any(ismember({'SET','INC'},type))
    error('makeLabel:invalidArguments','Must supply a valid type');
end
if ~isnumeric(value) && ~islogical(value)
    error('makeLabel:invalidArguments','Must supply a valid numeric or logical value');
end

switch(type)
    case 'SET'
        out.type = 'labelset';
    case 'INC'
        out.type = 'labelinc';
    otherwise
        error('unknown label type');
end
out.label=label;
out.value=value;
end
