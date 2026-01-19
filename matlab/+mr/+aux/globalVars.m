function out = globalVars(get_set,var,val)
%GLOBAL Auxilary function to get or set global variables

persistent SupportedLabels;

if strcmp(get_set,'get')    
    eval(['out=' var ';']);
elseif strcmp(get_set,'set')
    if nargin<3
        error('globalVars: missing the third argument');
    end    
    eval([var '=val;']);
else
    error('globalVars: the first parameter must be either get or set');
end

end

