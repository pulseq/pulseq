%!test %%% on Octave run with oruntests() %%%
%! testMakeRfShim
function tests = testMakeRfShim
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

%% Test basic shim vector
function test_basic_shim(testCase)
    shimVec = [1 0.8 0.9 1.1];
    rfShim = mr.makeRfShim(shimVec);
    testCase.verifyEqual(rfShim.type, 'rfShim');
    testCase.verifyEqual(rfShim.shimVector, shimVec(:), 'AbsTol', 1e-10);
end

%% Test row input becomes column
function test_row_to_column(testCase)
    shimVec = [1 0.5 0.7 0.3];
    rfShim = mr.makeRfShim(shimVec);
    testCase.verifyEqual(size(rfShim.shimVector, 2), 1, ...
        'shimVector should be a column vector');
end

%% Test column input stays column
function test_column_input(testCase)
    shimVec = [1; 0.5; 0.7];
    rfShim = mr.makeRfShim(shimVec);
    testCase.verifyEqual(rfShim.shimVector, shimVec, 'AbsTol', 1e-10);
end

%% Test complex shim values
function test_complex_shim(testCase)
    shimVec = [1+0.5i, 0.8-0.3i, 0.9+0.1i];
    rfShim = mr.makeRfShim(shimVec);
    testCase.verifyEqual(rfShim.shimVector, shimVec(:), 'AbsTol', 1e-10);
end

%% Test single element
function test_single_element(testCase)
    rfShim = mr.makeRfShim(1.5);
    testCase.verifyEqual(rfShim.shimVector, 1.5);
end
