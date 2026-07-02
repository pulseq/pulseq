%!test %%% on Octave run with oruntests() %%%
%! testMakeAdiabaticPulse
function tests = testMakeAdiabaticPulse
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

%% Test create hypsec pulse
function test_hypsec(testCase)
    testCase.assumeTrue(mr.aux.isSigPyAvailable(), ...
        'Python/sigpy not available, skipping adiabatic pulse tests');
    rf = mr.makeAdiabaticPulse('hypsec');
    testCase.verifyEqual(rf.type, 'rf');
    testCase.verifyTrue(length(rf.signal) > 0);
end

%% Test create wurst pulse
function test_wurst(testCase)
    testCase.assumeTrue(mr.aux.isSigPyAvailable(), ...
        'Python/sigpy not available, skipping adiabatic pulse tests');
    rf = mr.makeAdiabaticPulse('wurst','duration', 40e-3);
    testCase.verifyEqual(rf.type, 'rf');
    testCase.verifyTrue(length(rf.signal) > 0);
end

%% Test hypsec with custom duration
function test_hypsec_duration(testCase)
    testCase.assumeTrue(mr.aux.isSigPyAvailable(), ...
        'Python/sigpy not available, skipping adiabatic pulse tests');
    dur = 15e-3;
    rf = mr.makeAdiabaticPulse('hypsec', 'Duration', dur);
    testCase.verifyEqual(rf.shape_dur, dur, 'RelTol', 0.01);
end

%% Test with slice thickness returns gradient
function test_with_slice_thickness(testCase)
    testCase.assumeTrue(mr.aux.isSigPyAvailable(), ...
        'Python/sigpy not available, skipping adiabatic pulse tests');
    [rf, gz] = mr.makeAdiabaticPulse('hypsec', 'sliceThickness', 5e-3);
    testCase.verifyEqual(rf.type, 'rf');
    testCase.verifyEqual(gz.type, 'trap');
    testCase.verifyEqual(gz.channel, 'z');
end
