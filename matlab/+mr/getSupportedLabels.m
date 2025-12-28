function supported_labels = getSupportedLabels()
% auxilary function

supported_labels={ ...
    ... % data counters
    'SLC','SEG','REP','AVG','SET','ECO','PHS','LIN','PAR','ACQ', ... % these are copied to the corresponding MDH fields on Siemens
    ... % data flags
    'NAV','REV','SMS', ... % these are copied to the corresponding MDH fields on Siemens 
    'REF', 'IMA', ... % flags for parallel imaging
    'OFF', ... % Offline flag that labels the data, that should not be used for the online-reconstruction (on Siemens it negates the ONLINE MDH flag)
    'NOISE', ... % flag marking noise adjust scan, for parallel imaging acceleration
    ... % control flags/switches -- they are not affecting the data but rather the sequence itself
        'PMC', ... # for MoCo/PMC Pulseq version to recognize blocks that can be prospectively corrected for motion 
        'NOROT','NOPOS','NOSCL', ... # instruct the interpreter to ignore the position, rotation or scaling of the FOV specified on the UI 
        'ONCE', ... # a 3-state flag that instructs the interpreter to alter the sequence when executing multiple repeats as follows: blocks with ONCE==0 are executed on every repetition; ONCE==1: only the first repetition; ONCE==2: only the last repetition
        'TRID', ... # an integer ID of the TR (sequence segment) used by the GE interpreter (and some others) to optimize the execution on the scanner
	... % mgram: add labels for dictionary simulation
	'START', ... # mgram; indicates start position for the MRF .seq simulation
	'STOP' ... # mgram; indicates stop position for the MRF .seq simulation        
    };

end
