%!test %%% on Octave run with oruntests() %%%
%! testAuxStrstrip
function tests = testAuxStrstrip
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

%% Test leading and trailing spaces
function test_leading_trailing_spaces(testCase)
    out = mr.aux.strstrip('  hello  ');
    testCase.verifyEqual(out, 'hello');
end

%% Test no whitespace
function test_no_whitespace(testCase)
    out = mr.aux.strstrip('hello');
    testCase.verifyEqual(out, 'hello');
end

%% Test empty string
function test_empty(testCase)
    out = mr.aux.strstrip('');
    testCase.verifyTrue(isempty(out));
end

%% Test only spaces
function test_only_spaces(testCase)
    out = mr.aux.strstrip('   ');
    testCase.verifyTrue(isempty(out));
end

%% Test tabs and newlines
function test_tabs_newlines(testCase)
    in = sprintf('\t hello \n');
    out = mr.aux.strstrip(in);
    testCase.verifyEqual(out, 'hello');
end

%% Test preserves internal whitespace
function test_internal_whitespace(testCase)
    out = mr.aux.strstrip('  hello world  ');
    testCase.verifyEqual(out, 'hello world');
end
