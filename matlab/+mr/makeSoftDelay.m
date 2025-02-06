function sd = makeSoftDelay(varargin)
%makeSoftDelay Create a soft delay extension event.
%   softDelay=makeSoftDelay() Create a soft delay event, that can be used in 
%                         combionation with an empty (pure delay) block
%                         e.g. to adjust TE, TR or other delays. The soft delay 
%                         extension acts by rewriting the block duration
%                         based on the user input (to the interpreter)
%                         according to the equation dur=input/factor+offset. 
%                         Required parameters are 'numeric ID' and 'string
%                         hint'. Optional parameter 'factor' can be either
%                         positive and negative. Optional parameter 'offset' 
%                         given in seconds can also be either positive and
%                         negative. The 'hint' parameter is expected to be
%                         for identical 'numID'. 
%
%   See also  Sequence.addBlock() Sequence.applySoftDelay()

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeSoftDelay';
    
    addRequired(parser, 'numID', @isnumeric);
    addRequired(parser, 'hint', @(x) ischar(x)&&~isempty(x)); 
    addOptional(parser, 'offset', 0, @isnumeric); 
    addOptional(parser, 'factor', 1, @isnumeric); 
end

parse(parser, varargin{:});
opt = parser.Results;

if length(regexp(opt.hint, '(\s+)','split'))>1
    error('makeSoftDelay: parameter ''hint'' may not contain white space caharacters');
end

sd.type   = 'softDelay';
sd.num    = opt.numID;
sd.hint   = opt.hint;
sd.offset = opt.offset;
sd.factor = opt.factor;

end
