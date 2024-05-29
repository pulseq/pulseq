function supported_labels = getSupportedLabels()
% auxilary function

supported_labels={'SLC','SEG','REP','AVG','SET','ECO','PHS','LIN','PAR','NAV','REV','SMS', ... # these label values are copied to the corresponding MDH field on Siemens
    'REF', 'IMA', ... # For parallel imaging
    'NOISE', ... % noise adjust scan, for iPAT acceleration
        'PMC', ... # for MoCo/PMC Pulseq version to recognize blocks that can be prospectively corrected for motion 
        'NOROT','NOPOS','NOSCL', ... # instruct the interpreter to ignore the position, rotation or scaling of the FOV specified on the UI 
        'ONCE', ... # a 3-state flag that instructs the interpreter to alter the sequence when executing multiple repeats as follows: blocks with ONCE==0 are executed on every repetition; ONCE==1: only the first repetition; ONCE==2: only the last repetition
        'TRID' ... # an integer ID of the TR (sequence segment) used by the GE interpreter to optimize the execution on the scanner
    };

end
