function tests = testCheckTiming
    tests = functiontests(localfunctions);
end

function test_check_timing(testCase)
    import mr.*

    % Define system settings
    system = opts( ...
        'maxGrad', 28, ...
        'gradUnit', 'mT/m', ...
        'maxSlew', 200, ...
        'slewUnit', 'T/m/s', ...
        'rfRingdownTime', 20e-6, ...
        'rfDeadTime', 100e-6, ...
        'adcDeadTime', 10e-6 ...
    );

    % System with zero ringdown and dead times to introduce errors
    systemBroken = opts( ...
        'maxGrad', 28, ...
        'gradUnit', 'mT/m', ...
        'maxSlew', 200, ...
        'slewUnit', 'T/m/s', ...
        'rfRingdownTime', 0, ...
        'rfDeadTime', 0, ...
        'adcDeadTime', 0 ...
    );

    % Create a sequence
    seq = Sequence(system);

    % Add events with possible timing errors
    rf = makeSincPulse(1, 'duration', 1e-3, 'delay', system.rfDeadTime, 'system', system);
    seq.addBlock(rf); % Block 1: No error

    rf = makeSincPulse(1, 'duration', 1e-3, 'system', systemBroken);
    seq.addBlock(rf); % Block 2: RF_DEAD_TIME, RF_RINGDOWN_TIME, BLOCK_DURATION_MISMATCH

    adc = makeAdc(100, 'duration', 1e-3, 'delay', system.adcDeadTime, 'system', system);
    seq.addBlock(adc); % Block 3: No error

    adc = makeAdc(123, 'duration', 1e-3, 'delay', system.adcDeadTime, 'system', system);
    seq.addBlock(adc); % Block 4: RASTER

    adc = makeAdc(100, 'duration', 1e-3, 'system', systemBroken);
    seq.addBlock(adc); % Block 5: ADC_DEAD_TIME, POST_ADC_DEAD_TIME

    gx = makeTrapezoid('x', 'area', 1, 'duration', 1, 'system', system);
    seq.addBlock(gx); % Block 6: No error

    gx = makeTrapezoid('x', 'area', 1, 'duration', 1.00001e-3, 'system', system);
    seq.addBlock(gx); % Block 7: RASTER

    gx = makeTrapezoid('x', 'flatArea', 1, 'riseTime', 1e-6, 'flatTime', 1e-3, 'fallTime', 3e-6, 'system', system);
    seq.addBlock(gx); % Block 8: RASTER

    gx = makeTrapezoid('x', 'area', 1, 'duration', 1e-3, 'delay', -1e-5, 'system', system);
    seq.addBlock(gx); % Block 9: NEGATIVE_DELAY

    % Check timing errors
    [ok, errorReport] = seq.checkTiming();
    
    % Verify that certain blocks have no errors
    % testCase.verifyTrue(blocksNotInErrorReport(errorReport, [1, 3, 6]), 'No timing errors expected on blocks 1, 3, and 6');
    testCase.verifyFalse(existsInErrorReport(errorReport, 1, 'Block:1'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 2, 'Block:1'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 3, 'Block:1'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 4, 'Block:1'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 5, 'Block:1'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 6, 'Block:1'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 1, 'Block:3'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 2, 'Block:3'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 3, 'Block:3'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 4, 'Block:3'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 5, 'Block:3'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 6, 'Block:3'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 1, 'Block:6'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 2, 'Block:6'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 3, 'Block:6'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 4, 'Block:6'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 5, 'Block:6'));
    testCase.verifyFalse(existsInErrorReport(errorReport, 6, 'Block:6'));

    % Verify specific errors exist in the report
    testCase.verifyTrue(existsInErrorReport(errorReport, 1, 'Block:2'));
    testCase.verifyTrue(existsInErrorReport(errorReport, 1, 'RF dead time'));
    testCase.verifyTrue(existsInErrorReport(errorReport, 1, 'rfRingdownTime'));
    
    testCase.verifyTrue(existsInErrorReport(errorReport, 2, 'Block:4'));
    testCase.verifyTrue(existsInErrorReport(errorReport, 2, 'adcRasterTime'));

    testCase.verifyTrue(existsInErrorReport(errorReport, 3, 'Block:5'));
    testCase.verifyTrue(existsInErrorReport(errorReport, 3, 'adcDeadTime'));
    testCase.verifyTrue(existsInErrorReport(errorReport, 3, 'post-adc'));

    testCase.verifyTrue(existsInErrorReport(errorReport, 4, 'Block:7'));
    testCase.verifyTrue(existsInErrorReport(errorReport, 4, 'blockDurationRaster'));
    
    testCase.verifyTrue(existsInErrorReport(errorReport, 5, 'Block:8'));
    testCase.verifyTrue(existsInErrorReport(errorReport, 5, 'blockDurationRaster'));
    testCase.verifyTrue(existsInErrorReport(errorReport, 5, 'riseTime'));
    testCase.verifyTrue(existsInErrorReport(errorReport, 5, 'fallTime'));

    testCase.verifyTrue(existsInErrorReport(errorReport, 6, 'Block:9'));
    testCase.verifyTrue(existsInErrorReport(errorReport, 6, 'delay:-'));
    
end

% Helper Functions
function result = existsInErrorReport(errorReport, block, errorSubstring)
    % Check if a specific error message substring is present in the error report for a block
    if block > length(errorReport) || isempty(errorReport{block})
        result = false;
        return;
    end
    result = contains(errorReport{block}, errorSubstring, 'IgnoreCase', true);
end