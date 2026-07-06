%!test %%% on Octave run with oruntests() %%%
%! testRfOptionalDelayOutput
function tests = testRfOptionalDelayOutput
    try
        mr.opts();
    catch
        pulseqPath=fullfile(fileparts(mfilename),'..','matlab');
        addpath(genpath(pulseqPath));
    end
    if exist('functiontests')
        tests = functiontests(localfunctions);
    else
        lf=localfunctions();
        testCase=makeOctaveTestCase();
        for i=1:length(lf)
            f=lf{i};
            n=func2str(f);
            if length(n)>3 && strcmp(n(1:4),'test')
                f(testCase);
                fprintf('Test function %s completed successfully\n', n);
            end
        end
    end
end

function test_optional_delay_output_duration(testCase)
    sys = mr.opts('rfDeadTime', 80e-6, 'rfRingdownTime', 30e-6);
    ringdownTime = getRfRingdownTime(sys);

    signal = ones(1, 128);

    [rfArb, ~, ~, delayArb] = mr.makeArbitraryRf(signal, pi/6, sys, ...
        'bandwidth', 1500, 'sliceThickness', 5e-3, 'delay', 40e-6);
    verifyDelayDuration(testCase, rfArb, delayArb, ringdownTime);

    [rfBlock, delayBlock] = mr.makeBlockPulse(pi/2, sys, ...
        'duration', 1.2e-3, 'delay', 40e-6);
    verifyDelayDuration(testCase, rfBlock, delayBlock, ringdownTime);

    [rfGauss, ~, ~, delayGauss] = mr.makeGaussPulse(pi/2, sys, ...
        'duration', 1.4e-3, 'timeBwProduct', 3, 'sliceThickness', 5e-3, 'delay', 40e-6);
    verifyDelayDuration(testCase, rfGauss, delayGauss, ringdownTime);

    [rfSinc, ~, ~, delaySinc] = mr.makeSincPulse(pi/2, sys, ...
        'duration', 1.6e-3, 'timeBwProduct', 4, 'sliceThickness', 5e-3, 'delay', 40e-6);
    verifyDelayDuration(testCase, rfSinc, delaySinc, ringdownTime);

    if mr.aux.isSigPyAvailable()
        [rfAdiabatic, ~, ~, delayAdiabatic] = mr.makeAdiabaticPulse('hypsec', sys, ...
            'duration', 8e-3, 'sliceThickness', 5e-3, 'delay', 40e-6);
        verifyDelayDuration(testCase, rfAdiabatic, delayAdiabatic, ringdownTime);

        [rfSLR, ~, ~, delaySLR] = mr.makeSLRpulse(pi/2, sys, ...
            'duration', 2e-3, 'timeBwProduct', 4, 'sliceThickness', 5e-3, 'delay', 40e-6);
        verifyDelayDuration(testCase, rfSLR, delaySLR, ringdownTime);
    else
        testCase.log('Skipping makeAdiabaticPulse/makeSLRpulse delay checks: sigpy not available\n');
    end
end

function verifyDelayDuration(testCase, rf, delayObj, ringdownTime)
    expectedDuration = rf.delay + rf.shape_dur + ringdownTime;
    testCase.verifyEqual(delayObj.type, 'delay');
    testCase.verifyEqual(delayObj.delay, expectedDuration, 'AbsTol', 1e-12);
end

function ringdownTime = getRfRingdownTime(sys)
    if isfield(sys, 'rfRingDownTime')
        ringdownTime = sys.rfRingDownTime;
    else
        ringdownTime = sys.rfRingdownTime;
    end
end
