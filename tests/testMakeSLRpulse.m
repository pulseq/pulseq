%!test %%% on Octave run with oruntests() %%%
%! testMakeSLRpulse
function tests = testMakeSLRpulse
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

%% Test basic SLR pulse creation
function test_basic_slr(testCase)
    testCase.assumeTrue(mr.aux.isSigPyAvailable(), ...
        'Python/sigpy not available, skipping SLR pulse tests');
    rf = mr.makeSLRpulse(pi/2, 'Duration', 3e-3, 'timeBwProduct', 4);
    testCase.verifyEqual(rf.type, 'rf');
    testCase.verifyTrue(length(rf.signal) > 0);
end

%% Test SLR with slice thickness returns gradient
function test_slr_with_slice(testCase)
    testCase.assumeTrue(mr.aux.isSigPyAvailable(), ...
        'Python/sigpy not available, skipping SLR pulse tests');
    [rf, gz] = mr.makeSLRpulse(pi/2, 'Duration', 3e-3, ...
        'timeBwProduct', 4, 'sliceThickness', 5e-3);
    testCase.verifyEqual(rf.type, 'rf');
    testCase.verifyEqual(gz.type, 'trap');
end

%% Test SLR excitation pulse
function test_slr_excitation(testCase)
    testCase.assumeTrue(mr.aux.isSigPyAvailable(), ...
        'Python/sigpy not available, skipping SLR pulse tests');
    rf = mr.makeSLRpulse(pi/2, 'Duration', 3e-3, ...
        'timeBwProduct', 4, 'use', 'excitation');
    testCase.verifyEqual(rf.use, 'excitation');
end
