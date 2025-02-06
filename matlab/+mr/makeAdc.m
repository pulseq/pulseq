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
%   given in sys. For example, adcDeadTime can be taken into account by 
%   making sure that both the before the ADC is as least as long as 
%   adcDeadTime and another period of adcDeadTime is added after sampling 
%   is finished. 
%
%   See also  Sequence.addBlock

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeAdc';
    
    addRequired(parser,'numSamples',@(x)(isnumeric(x) && (fix(x)-x)==0));
    addOptional(parser,'system',[],@isstruct);
    addParamValue(parser,'dwell',0,@isnumeric);
    addParamValue(parser,'duration',0,@isnumeric);
    addParamValue(parser,'delay',0,@isnumeric);
    addParamValue(parser,'ppmOffset', 0, @isnumeric);
    addParamValue(parser,'freqOffset',0,@isnumeric);
    addParamValue(parser,'phaseOffset',0,@isnumeric);
    addParamValue(parser,'phaseModulation',[],@isnumeric);
end

parse(parser,num,varargin{:});
opt = parser.Results;

if isempty(opt.system)
    system=mr.opts();
else
    system=opt.system;
end

adc.type = 'adc';
adc.numSamples = num;
adc.dwell = opt.dwell;
adc.delay = opt.delay;
adc.ppmOffset = opt.ppmOffset;
adc.freqOffset = opt.freqOffset;
adc.phaseOffset = opt.phaseOffset;
adc.deadTime = system.adcDeadTime;

if (opt.dwell==0 && opt.duration==0) || (opt.dwell>0 && opt.duration>0)
    error('Either dwell or duration must be defined');
end

if ~isempty(opt.phaseModulation)
    if length(opt.phaseModulation)~=num
        error('ADC Phase modulation vector must have the same length as the number of samples');
    end
    adc.phaseModulation=opt.phaseModulation;
else
    adc.phaseModulation=[];
end

if opt.duration > 0
    adc.dwell = opt.duration/opt.numSamples;
end
if opt.dwell > 0
    adc.duration = opt.dwell*opt.numSamples;
end
if adc.deadTime > adc.delay
    adc.delay = adc.deadTime; % adcDeadTime is added before the actual sampling (and also second time after the sampling period)
end

end