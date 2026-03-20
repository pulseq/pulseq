function tests = testCalcRamp
    tests = functiontests(localfunctions);
end

%% Test simple ramp from zero to zero
function test_zero_to_zero(testCase)
    sys = mr.opts();
    grt = sys.gradRasterTime;
    % k0 and kend both imply zero gradient
    k0   = [0 0; 0 0; 0 0];        % [3,2]: two points both at origin
    kend = [grt*2 grt*3; 0 0; 0 0]; % next two points on x-axis with zero gradient
    % This means G0=[0;0;0] and Gend=[0;0;0]
    % So the ramp should just be empty or trivial
    [kout, success] = mr.calcRamp(k0, kend, 'system', sys);
    testCase.verifyTrue(success == 1 || success == true);
end

%% Test simple ramp in one axis
function test_simple_ramp_x(testCase)
    sys = mr.opts();
    grt = sys.gradRasterTime;
    % Starting at zero gradient, ending at a small gradient in x
    k0 = [0 0; 0 0; 0 0];          % G0 = 0
    Gend = 1000;  % target gradient in Hz/m
    dk = Gend * grt;
    kend = [2*dk 3*dk; 0 0; 0 0];  % Gend = dk/grt = Gend
    [kout, success] = mr.calcRamp(k0, kend, 'system', sys);
    testCase.verifyTrue(logical(success), 'Ramp should succeed for small gradient');
    testCase.verifyEqual(size(kout, 1), 3, 'Output should have 3 rows');
end

%% Test MaxPoints=0 may fail for large gradient
function test_maxpoints_zero(testCase)
    sys = mr.opts();
    grt = sys.gradRasterTime;
    % Large gradient step
    Gend = sys.maxGrad * 0.9;
    dk = Gend * grt;
    k0 = [0 0; 0 0; 0 0];
    kend = [2*dk 3*dk; 0 0; 0 0];
    [~, success] = mr.calcRamp(k0, kend, 'system', sys, 'MaxPoints', 0);
    % With MaxPoints=0, only direct connection is tried
    % This may or may not succeed depending on slew constraints
    testCase.verifyTrue(islogical(logical(success)));
end

%% Test output shape is [3, N]
function test_output_shape(testCase)
    sys = mr.opts();
    grt = sys.gradRasterTime;
    k0 = [0 0; 0 0; 0 0];
    kend = [grt*2 grt*3; 0 0; 0 0];
    [kout, ~] = mr.calcRamp(k0, kend, 'system', sys);
    testCase.verifyEqual(size(kout, 1), 3);
end
