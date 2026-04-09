function adc=makeAdc(num,varargin)
%makeAdc Create an ADC readout event.
%
%   PURPOSE
%     Build an ADC (analog-to-digital converter) sampling event struct.
%     The returned struct describes when sampling occurs, how many
%     samples are acquired, and at what dwell rate, and is consumed by
%     mr.Sequence/addBlock to add the readout to a sequence.
%
%   SIGNATURES
%     adc = mr.makeAdc(numSamples, 'Duration', d)
%     adc = mr.makeAdc(numSamples, 'Dwell', dt)
%     adc = mr.makeAdc(numSamples, system, ...)              % system as 2nd positional arg
%     adc = mr.makeAdc(numSamples, ..., 'system', system)    % system as name/value
%     adc = mr.makeAdc(numSamples, ..., 'Delay', d)
%
%     Exactly one of 'Duration' or 'Dwell' must be supplied.
%     If system is omitted, mr.opts() is used.
%     Parameter names are case-insensitive.
%
%   INPUTS
%     numSamples         integer  number of ADC samples (must be a whole number)  required
%     system             struct   from mr.opts; defaults to mr.opts() if omitted  optional
%     'Duration'         double   total acquisition duration, seconds. Dwell is
%                                 computed as Duration/numSamples.                name/value
%     'Dwell'            double   per-sample dwell time, seconds. Total acquisition
%                                 time is Dwell*numSamples.                       name/value
%     'Delay'            double   delay before sampling starts, seconds, default 0.
%                                 Silently bumped to system.adcDeadTime if smaller
%                                 (see NOTES).                                    name/value
%     'freqOffset'       double   demodulation frequency offset, Hz, default 0    name/value
%     'phaseOffset'      double   demodulation phase offset, radians, default 0   name/value
%     'freqPPM'          double   frequency offset in PPM (relative to system B0
%                                 and gamma), default 0                           name/value
%     'phasePPM'         double   phase offset in PPM, default 0                  name/value
%     'phaseModulation'  vector   per-sample phase modulation, length must equal
%                                 numSamples. Default: [] (no modulation).        name/value
%
%   OUTPUT
%     adc  struct with fields (in order returned by fieldnames):
%       .type             char,    always 'adc'
%       .numSamples       integer, number of samples (= input numSamples)
%       .delay            double,  delay before sampling, seconds
%                                  (>= system.adcDeadTime; see NOTES)
%       .freqOffset       double,  Hz
%       .phaseOffset      double,  radians
%       .freqPPM          double,  PPM
%       .phasePPM         double,  PPM
%       .deadTime         double,  copy of system.adcDeadTime at construction, seconds
%       .phaseModulation  double vector or []
%       .dwell            double,  sample dwell time, seconds
%                                  (= Duration/numSamples if Duration was given,
%                                   else the supplied Dwell)
%
%   ERRORS
%     - 'Either dwell or duration must be defined': both are 0, or both > 0.
%       Exactly one must be specified.
%     - 'ADC Phase modulation vector must have the same length as the number
%       of samples': length(phaseModulation) ~= numSamples.
%     - MATLAB:InputParser:* errors: numSamples is not an integer, or other
%       validator failures (non-numeric where numeric is expected, etc.).
%
%   NOTES
%     - The .delay field is silently increased to system.adcDeadTime if the
%       requested delay is smaller. If you set 'Delay' to 0 with a nonzero
%       system.adcDeadTime, the returned adc.delay will equal system.adcDeadTime,
%       not 0. Account for this when computing block timing.
%     - system.adcDeadTime is captured into adc.deadTime at construction time;
%       changing system later does not retroactively update existing adc events.
%     - Caches an inputParser in a persistent variable for performance;
%       no other global state.
%
%   EXAMPLE
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s', ...
%                   'adcDeadTime', 20e-6);
%     Nx = 256;  fov = 256e-3;  deltak = 1/fov;
%     gx  = mr.makeTrapezoid('x', sys, 'FlatArea', Nx*deltak, 'FlatTime', 6.4e-3);
%     % ADC sampled across the readout flat top
%     adc = mr.makeAdc(Nx, sys, 'Duration', gx.flatTime, 'Delay', gx.riseTime);
%
%   SEE ALSO
%     mr.opts, mr.makeTrapezoid, mr.Sequence/addBlock, mr.calcDuration

persistent parser
if isempty(parser)
    parser = mr.aux.InputParserCompat;
    parser.FunctionName = 'makeAdc';
    
    addRequired(parser,'numSamples',@(x)(isnumeric(x) && (fix(x)-x)==0));
    addOptional(parser,'system',[],@isstruct);
    addParamValue(parser,'dwell',0,@isnumeric);
    addParamValue(parser,'duration',0,@isnumeric);
    addParamValue(parser,'delay',0,@isnumeric);
    addParamValue(parser,'freqOffset',0,@isnumeric);
    addParamValue(parser,'phaseOffset',0,@isnumeric);
    addParamValue(parser,'freqPPM', 0, @isnumeric);
    addParamValue(parser,'phasePPM', 0, @isnumeric);
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
adc.delay = opt.delay;
adc.freqOffset = opt.freqOffset;
adc.phaseOffset = opt.phaseOffset;
adc.freqPPM = opt.freqPPM;
adc.phasePPM = opt.phasePPM;
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
else
    adc.dwell = opt.dwell;
end
%if opt.dwell > 0
%    adc.duration = opt.dwell*opt.numSamples;
%end
if adc.deadTime > adc.delay
    adc.delay = adc.deadTime; % adcDeadTime is added before the actual sampling (and also second time after the sampling period)
end

end