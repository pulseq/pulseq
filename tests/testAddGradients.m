%!test %%% on Octave run with oruntests() %%%
%! testAddGradients
function tests = testAddGradients
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

%% Test sum of two identical trapezoids doubles amplitude
function test_sum_identical_traps(testCase)
    g = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    gsum = mr.addGradients({g, g});
    testCase.verifyEqual(gsum.amplitude, 2*g.amplitude, 'AbsTol', 1e-6);
    testCase.verifyEqual(gsum.area, 2*g.area, 'AbsTol', 1e-3);
    testCase.verifyEqual(gsum.flatArea, 2*g.flatArea, 'AbsTol', 1e-3);
end

%% Test opposing trapezoids cancel
function test_opposing_traps(testCase)
    g1 = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    g2 = mr.makeTrapezoid('x', 'Area', -1000, 'Duration', 1e-3);
    gsum = mr.addGradients({g1, g2});
    testCase.verifyEqual(gsum.area, 0, 'AbsTol', 1);
end

%% Test non-cell input error
function test_non_cell_error(testCase)
    g = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    testCase.verifyError(@() mr.addGradients(g),'');
end

%% Test single gradient error
function test_single_gradient_error(testCase)
    g = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    testCase.verifyError(@() mr.addGradients({g}),'');
end

%% Test different channels error
function test_different_channels_error(testCase)
    gx = mr.makeTrapezoid('x', 'Area', 1000, 'Duration', 1e-3);
    gy = mr.makeTrapezoid('y', 'Area', 1000, 'Duration', 1e-3);
      testCase.verifyError(@() mr.addGradients({gx, gy}),'');
end

%% Test trap + extended trap
function test_trap_plus_extended(testCase)
    sys = mr.opts();
    g_trap = mr.makeTrapezoid('x', 'Amplitude', 1e5, 'Duration', 1e-3);
    g_ext = mr.makeExtendedTrapezoid('x', ...
        'Times', [0 2e-4 8e-4 1e-3], ...
        'Amplitudes', [0 50000 50000 0]);
    gsum = mr.addGradients({g_trap, g_ext}, sys);
    % Result should be an extended trapezoid
    testCase.verifyTrue(strcmp(gsum.type, 'grad'));
    % Duration should be the max of both
    dur_sum = mr.calcDuration(gsum);
    dur_max = max(mr.calcDuration(g_trap), mr.calcDuration(g_ext));
    testCase.verifyEqual(dur_sum, dur_max, 'AbsTol', sys.gradRasterTime);
end

%% Test three trapezoids
function test_three_traps(testCase)
    sys = mr.opts();
    g1 = mr.makeTrapezoid('y', 'Area', 1000, 'system', sys); 
    g2 = mr.makeTrapezoid('y', 'Area', 700, 'system', sys);
    g3 = mr.makeTrapezoid('y', 'Area', 500, 'system', sys);
    g4 = mr.makeTrapezoid('y', 'Area', 1000, 'system', sys, 'maxGrad', sys.maxGrad/3); % derate the amplitude to avoid exceeding slew rate limits due to simultaneous slewing
    g5 = mr.makeTrapezoid('y', 'Area', 700, 'system', sys, 'maxGrad', sys.maxGrad/3);
    g6 = mr.makeTrapezoid('y', 'Area', 500, 'system', sys, 'maxGrad', sys.maxGrad/3);
    g7 = mr.makeTrapezoid('y', 'Area', 1000, 'system', sys, 'maxSlew', sys.maxSlew/3); % derate the slew to avoid exceeding slew rate limits due to simultaneous slewing
    g8 = mr.makeTrapezoid('y', 'Area', 700, 'system', sys, 'maxSlew', sys.maxSlew/3);
    g9 = mr.makeTrapezoid('y', 'Area', 500, 'system', sys, 'maxSlew', sys.maxSlew/3);
    ga = mr.makeTrapezoid('y', 'Area', 1000, 'system', sys, 'maxGrad', sys.maxGrad/3, 'maxSlew', sys.maxSlew/3); % derate both the max grad slew to avoid exceeding slew rate limits due to simultaneous slewing
    gb = mr.makeTrapezoid('y', 'Area', 700, 'system', sys, 'maxGrad', sys.maxGrad/3, 'maxSlew', sys.maxSlew/3);
    gc = mr.makeTrapezoid('y', 'Area', 500, 'system', sys, 'maxGrad', sys.maxGrad/3, 'maxSlew', sys.maxSlew/3);
    testCase.verifyError(@() mr.addGradients({g1, g2, g3}),''); % must fail due to exceeding both amplitude and slew rate limits
    testCase.verifyError(@() mr.addGradients({g4, g5, g6}),''); % must fail due to exceeding slew rate limits
    testCase.verifyError(@() mr.addGradients({g7, g8, g9}),''); % must fail due to exceeding amplitude limits
    gsum = mr.addGradients({ga, gb, gc});
    testCase.verifyEqual(gsum.area, ga.area + gb.area + gc.area, 'AbsTol', 1);
end
