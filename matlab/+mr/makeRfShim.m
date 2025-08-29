function rfShim = makeRfShim(shimVec)
%makeRfShim Create a RF shimming event.
%   rfShim=makeRfShim(shimVec) Create an event describing RF shimmming in 
%                         the current block. It is only valid if the block 
%                         contains an RF pulse. 
%
%   See also  Sequence.addBlock

rfShim=struct('type','rfShim','shimVector',shimVec);

end
