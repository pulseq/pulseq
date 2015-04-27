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
%   See also  Sequence.addBlock

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeAdc';
    
    addRequired(parser,'numSamples',@(x)(isnumeric(x) && (fix(x)-x)==0));
    addParamValue(parser,'dwell',0,@isnumeric);
    addParamValue(parser,'duration',0,@isnumeric);
    addParamValue(parser,'delay',0,@isnumeric);
    addParamValue(parser,'freqOffset',0,@isnumeric);
    addParamValue(parser,'phaseOffset',0,@isnumeric);
end

parse(parser,num,varargin{:});

adc = parser.Results;
adc.type = 'adc';

if (adc.dwell==0 && adc.duration==0) || (adc.dwell>0 && adc.duration>0)
    error('Either dwell or duration must be defined');
end

if adc.duration > 0
    adc.dwell = adc.duration/adc.numSamples;
end
if adc.dwell > 0
    adc.duration = adc.dwell*adc.numSamples;
end
end