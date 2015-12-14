function adc=makeAdc(num,varargin)
%makeAdc Create a ADC readout event.
%   adc=makeAdc(n, 'Dwell', dwell) Create ADC with n samples
%   with given dwell time.
%
%   adc=makeAdc(n, 'Duration', dur) Create ADC with n
%   samples and specified total duration.
%
%   adc=makeAdc(..., 'Delay', d) Create ADC with initial delay.
%
%   adc=makeAdc(n, sys, ...) Create ADC considering system properties
%   given in sys. For example, a dead time after sampling can be added to
%   the duration. 
%
%   See also  Sequence.addBlock

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeAdc';
    
    addRequired(parser,'numSamples',@(x)(isnumeric(x) && (fix(x)-x)==0));
    addOptional(parser,'system',mr.opts(),@isstruct);
    addParamValue(parser,'dwell',0,@isnumeric);
    addParamValue(parser,'duration',0,@isnumeric);
    addParamValue(parser,'delay',0,@isnumeric);
    addParamValue(parser,'freqOffset',0,@isnumeric);
    addParamValue(parser,'phaseOffset',0,@isnumeric);
end

parse(parser,num,varargin{:});

opt = parser.Results;
adc.type = 'adc';
adc.numSamples = num;
adc.dwell = opt.dwell;
adc.delay = opt.delay;
adc.freqOffset = opt.freqOffset;
adc.phaseOffset = opt.phaseOffset;
adc.deadTime = opt.system.adcDeadTime;

if (opt.dwell==0 && opt.duration==0) || (opt.dwell>0 && opt.duration>0)
    error('Either dwell or duration must be defined');
end

if opt.duration > 0
    adc.dwell = opt.duration/opt.numSamples;
end
if opt.dwell > 0
    adc.duration = opt.dwell*opt.numSamples;
end

end