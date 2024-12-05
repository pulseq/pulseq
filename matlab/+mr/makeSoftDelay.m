function sd = makeSoftDelay(varargin)
%makeSoftDelay Create a soft delay extension event.
%   softDelay=makeSoftDelay() Create a soft delay event, that can be used in 
%                         combionation with an empty (pure delay) block
%                         e.g. to adjust TE, TR or other delays. Required
%                         parameters are 'numeric ID' and 'string hint'.
%                         Optional parameter 'offset' given in seconds,
%                         which can be both positive and negative, will be
%                         added to the user-provided delay. The 'hint'
%                         parameter is not allowed to contain white spaces. 
%
%   See also  Sequence.addBlock

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeSoftDelay';
    
    addRequired(parser, 'number', @isnumeric);
    addRequired(parser, 'hint', @(x) ischar(x)&&~isempty(x)); 
    addOptional(parser, 'offset', 0, @isnumeric); 
end

parse(parser, varargin{:});
opt = parser.Results;

if length(regexp(opt.hint, '(\s+)','split'))>1
    error('makeSoftDelay: parameter ''hint'' may not contain white space caharacters');
end

sd.type   = 'softDelay';
sd.num    = opt.number;
sd.hint   = opt.hint;
sd.offset = opt.offset;

end
