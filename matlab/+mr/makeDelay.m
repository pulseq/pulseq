function del = makeDelay(delay)
%makeDelay Create a delay event.
%   delay=makeDelay(del) Create delay event with given delay del.
%
%   See also  Sequence.addBlock

if nargin<1
    error('makeDelay:invalidArguments','Must supply a delay');
end
assert(isfinite(delay) & delay>0,'makeDelay:invalidDelay',...
    'Delay (%.2f ms) is invalid',delay*1e3);
del.type = 'delay';
del.delay = delay;
end
