function tests = testAuxStrstrip
    tests = functiontests(localfunctions);
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
