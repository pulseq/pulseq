%!test %%% on Octave run with oruntests() %%%
%! testAuxVersion
function tests = testAuxVersion
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

%% Test default version string format
function test_default_format(testCase)
    [version_major, version_minor, version_revision, version_combined] = mr.aux.version();
    testCase.verifyTrue(version_major>=1);
    testCase.verifyTrue(version_minor>=5);
    testCase.verifyTrue(version_revision>=0);
    testCase.verifyTrue(version_combined>=1005001);
end

%% Test 'output' type returns struct with major, minor, revision
function test_output_struct(testCase)
    [version_major, version_minor, version_revision, version_combined] = mr.aux.version('output');
    testCase.verifyTrue(version_major>=1);
    testCase.verifyTrue(version_minor>=5);
    testCase.verifyTrue(version_revision>=0);
    testCase.verifyTrue(version_combined>=1005001);
end

%% Test combined version number
function test_combined(testCase)
    [version_major, version_minor, version_revision, version_combined] = mr.aux.version('output');
    expected_combined = version_major * 1000000 + version_minor * 1000 + version_revision;
    testCase.verifyEqual(version_combined, expected_combined);
end

