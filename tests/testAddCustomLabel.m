%!test %%% on Octave run with oruntests() %%%
%! testAddCustomLabel
function tests = testAddCustomLabel
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

%% Setup: reset global state
function setup(~)
    mr.aux.globalVars('set', 'SupportedLabels', []);
end

%% Teardown: reset global state
function teardown(~)
    mr.aux.globalVars('set', 'SupportedLabels', []);
end

%% Test add valid new label
function test_add_valid_label(testCase)
    mr.addCustomLabel('MYCUSTOM');
    labels = mr.getSupportedLabels();
    testCase.verifyTrue(any(strcmp(labels, 'MYCUSTOM')), ...
        'Custom label should appear in getSupportedLabels()');
end

% %% Test add duplicate label produces warning
% function test_add_duplicate_warning(testCase)
%     mr.addCustomLabel('TESTDUP');
%     % Second call with same label should warn
%     testCase.verifyWarning(@() mr.addCustomLabel('TESTDUP'));
% end

%% Test non-char input throws error
function test_non_char_error(testCase)
    %try
    %    mr.addCustomLabel(123); % should throw an error
    %    res=false;
    %catch
    %    res=true;
    %end
    %testCase.verifyTrue(res);
    testCase.verifyError(@() mr.addCustomLabel(123),'');
    %testCase.verifyError(@() mr.addCustomLabel(123), ?MException);
end

%% Test multiple custom labels
function test_multiple_labels(testCase)
    mr.addCustomLabel('CUSTOM1');
    mr.addCustomLabel('CUSTOM2');
    labels = mr.getSupportedLabels();
    testCase.verifyTrue(any(strcmp(labels, 'CUSTOM1')));
    testCase.verifyTrue(any(strcmp(labels, 'CUSTOM2')));
end

%% Test custom label works with makeLabel
function test_custom_label_with_makeLabel(testCase)
    mr.addCustomLabel('MYTEST');
    lbl = mr.makeLabel('SET', 'MYTEST', 42);
    testCase.verifyEqual(lbl.label, 'MYTEST');
    testCase.verifyEqual(lbl.value, 42);
end
