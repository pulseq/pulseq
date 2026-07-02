%!test %%% on Octave run with oruntests() %%%
%! testMakeSoftDelay
function tests = testMakeSoftDelay
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

%% Test valid creation
function test_valid_creation(testCase)
    sd = mr.makeSoftDelay(1, 'TE');
    testCase.verifyEqual(sd.type, 'softDelay');
    testCase.verifyEqual(sd.num, 1);
    testCase.verifyEqual(sd.hint, 'TE');
    testCase.verifyEqual(sd.offset, 0);
    testCase.verifyEqual(sd.factor, 1);
end

%% Test with custom offset and factor
function test_custom_offset_factor(testCase)
    sd = mr.makeSoftDelay(2, 'TR', -0.005, 2);
    testCase.verifyEqual(sd.num, 2);
    testCase.verifyEqual(sd.hint, 'TR');
    testCase.verifyEqual(sd.offset, -0.005);
    testCase.verifyEqual(sd.factor, 2);
end

%% Test hint with whitespace throws error
function test_whitespace_hint_error(testCase)
    testCase.verifyError(@() mr.makeSoftDelay(1, 'my delay'),'');
end

%% Test different numIDs
function test_different_numids(testCase)
    sd1 = mr.makeSoftDelay(1, 'TE');
    sd2 = mr.makeSoftDelay(99, 'TR');
    testCase.verifyEqual(sd1.num, 1);
    testCase.verifyEqual(sd2.num, 99);
end

%% Test negative factor
function test_negative_factor(testCase)
    sd = mr.makeSoftDelay(1, 'delay', 0, -1);
    testCase.verifyEqual(sd.factor, -1);
end
