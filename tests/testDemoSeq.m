% simple test that runs all demo sequences
function tests = testDemoSeq
    tests = functiontests(localfunctions);
end

function test_run_all_demo_seqs(testCase)
    % excluded sequences (too slow or known not to run)
    excludedSeqs={'writeMPRAGE_4ge', 'just_some_nonexistant_sequence_to_test_multiple_excludes'};
    % prepare environment and enumerate files
    currDir  = fileparts(which('testDemoSeq'));
    % test whether the main toolbox is in the path
    try
        sys=mr.opts();
    catch
        addpath(fullfile(currDir,'../matlab'));
    end
    demoDir=fullfile(currDir,'../matlab/demoSeq');
    d = dir(fullfile(demoDir, '*.m'));
    demoSeqs={};
    for k = 1:numel(d)
        [~, name, ~] = fileparts(d(k).name);
        safeName = matlab.lang.makeValidName(name);
        if any(strcmp(safeName,excludedSeqs))
            continue;
        end
        demoSeqs{end+1}=safeName;
    end
    addpath(fullfile(currDir,'../matlab/demoUnsorted')); % for fidBeep
    addpath(demoDir);
    % now run all sequences one after another
    succeededSeqs={};
    failedSeqs={};
    for currSeq=demoSeqs
        try
            fprintf('***** running sequence: %s *****\n', currSeq{1});
            eval(currSeq{1});
            succeededSeqs{end+1}=currSeq{1};
            close all;
        catch ME
            %testCase.verifyFail(...
            %    sprintf('Script %s threw: %s', s, ME.message));
            fprintf('failed sequence: %s\n', currSeq{1});
            fprintf('error message: %s\n', ME.message);
            failedSeqs{end+1}=['failed sequence: ''' currSeq{1} ''', error message: ' ME.message ' '];
        end
    end
    testCase.log(1, ['Succeded sequences:' sprintf(' %d', succeededSeqs{:})]);
    if isempty(failedSeqs)
        testCase.verifyTrue(true);
    else
        testCase.verifyFail(['Failed sequences: ' failedSeqs{:}]);
    end
end
