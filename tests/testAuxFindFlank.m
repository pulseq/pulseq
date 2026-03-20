function tests = testAuxFindFlank
    tests = functiontests(localfunctions);
end

%% Test symmetric Gaussian flank detection
function test_gaussian_flank(testCase)
    x = linspace(-3, 3, 601);
    y = exp(-x.^2 / (2*0.5^2));
    % flank at cutoff=0.5 should be on the left side
    xf = mr.aux.findFlank(x, y, 0.5);
    testCase.verifyGreaterThan(xf, x(1));
    testCase.verifyLessThan(xf, 0); % left flank is negative x
end

%% Test step function
function test_step_function(testCase)
    x = 1:10;
    y = [0 0 0 0 0 1 1 1 1 1];
    xf = mr.aux.findFlank(x, y, 0.5);
    testCase.verifyGreaterThanOrEqual(xf, 5);
    testCase.verifyLessThanOrEqual(xf, 6);
end

%% Test cutoff zero returns first point
function test_cutoff_zero(testCase)
    x = 1:5;
    y = [0.1 0.5 1 0.5 0.1];
    xf = mr.aux.findFlank(x, y, 0);
    testCase.verifyEqual(xf, x(1), 'AbsTol', 0.01);
end

%% Test plateau signal
function test_plateau(testCase)
    x = 1:5;
    y = [1 1 1 1 1];
    xf = mr.aux.findFlank(x, y, 0.5);
    testCase.verifyEqual(xf, x(1), 'AbsTol', 0.01);
end
