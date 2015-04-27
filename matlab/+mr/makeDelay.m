function del = makeDelay(delay)
%makeDelay Create a delay event.
%   delay=makeDelay(del) Create delay event with given delay del.
%
%   See also  Sequence.addBlock
if nargin<1
    error('makeDelay:invalidArguments','Must supply a delay');
end
del.type = 'delay';
del.delay = delay;
end
