% simple test that runs all demo sequences
function tests = testDemoSeq
  if exist('functiontests')
    tests = functiontests(localfunctions);
  else
    lf=localfunctions();
    for i=1:length(lf)
      f=lf{i};
      n=func2str(f);
      if length(n)>3 && strcmp(n(1:4),'test')
        f();
      end
    end
  end
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
    
    if mr.aux.isOctave
      excludedSeqs{end+1}='writeZTE_Petra'; % it works but takes too long to generate (like one day?)
      excludedSeqs{end+1}='writeZTE_Petra_sodium'; % even longer...
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
    finishedSeconds=[];
    for currSeq=demoSeqs
        try
            fprintf('***** running sequence: %s *****\n', currSeq{1});
            tstart=tic;
            eval(currSeq{1});
            succeededSeqs{end+1}=currSeq{1};
        catch ME
            %testCase.verifyFail(...
            %    sprintf('Script %s threw: %s', s, ME.message));
            fprintf('failed sequence: %s\n', currSeq{1});
            fprintf('error message: %s\n', ME.message);
            failedSeqs{end+1}=['failed sequence: ''' currSeq{1} ''', error message: ' ME.message ' '];
        end
        finishedSeconds(end+1)=toc(tstart);
        fprintf('===== sequence %s finished after %.3g seconds =====\n', currSeq{1}, finishedSeconds(end));
        close all;
        mr.opts('resetDefault',true);
        mr.aux.globalVars('reset','SupportedLabels');
    end

    if isempty(failedSeqs)
      for i=1:length(succeededSeqs)
        succeededSeqs{i}=[succeededSeqs{i} sprintf('(%.3gs)',finishedSeconds(i))];
      end
    end

    if mr.aux.isOctave()
      fprintf('%s\n',['Succeded sequences:' sprintf(' %s', succeededSeqs{:})]);
      if ~isempty(failedSeqs)
        fprintf('%s\n',['Failed sequences: ' failedSeqs{:}]);
      end
    else
      if isempty(failedSeqs)
        testCase.log(1, ['Succeded sequences:' sprintf(' %s', [succeededSeqs{:}])]);
        testCase.verifyTrue(true);
      else
        testCase.log(1, ['Succeded sequences:' sprintf(' %s', succeededSeqs{:})]);
        testCase.verifyFail(['Failed sequences: ' failedSeqs{:}]);
      end
    end
end
