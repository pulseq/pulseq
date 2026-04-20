% simple test that runs legacy integration test scripts.
%
%   Enumerates *.out files in tests/legacy/approved/ and creates one test
%   point per file.  Each test runs the corresponding script, then compares
%   the freshly generated .out file against the approved reference, ignoring
%   line-ending differences.
%
%   Run with:   runtests('testLegacy')

function tests = testLegacy
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

function test_enumerate_and_run_legacy_tests(testCase)
         % ---- local helper -----
        function result = compareTextFiles(file1, file2)
        %compareTextFiles  Compare two text files ignoring line-ending differences.
            fid1 = fopen(file1, 'rb');
            assert(fid1 ~= -1, 'Cannot open file: %s', file1);
            txt1 = fread(fid1, '*char')';
            fclose(fid1);

            fid2 = fopen(file2, 'rb');
            assert(fid2 ~= -1, 'Cannot open file: %s', file2);
            txt2 = fread(fid2, '*char')';
            fclose(fid2);

            CR = char(13);  LF = char(10);
            txt1 = strrep(txt1, [CR LF], LF);
            txt1 = strrep(txt1, CR, LF);
            txt2 = strrep(txt2, [CR LF], LF);
            txt2 = strrep(txt2, CR, LF);

            result = strcmp(txt1, txt2);
        end

    % enumerateFiles()
    currDir  = fileparts(which('testLegacy'));
    % test whether the main toolbox is in the path
    try
        sys=mr.opts();
    catch
        addpath(fullfile(currDir,'../matlab'));
    end
    approvedDir = fullfile(currDir, 'legacy', 'approved');
    d = dir(fullfile(approvedDir, '*.matlab.out'));
    files = {};
    for k = 1:numel(d)
        [~, name, ~] = fileparts(d(k).name);
        if strcmp(name(end-6:end),'.matlab'), name=name(1:end-7); end
        files{end+1} = name;
    end

    %verifyLegacyFile(testCase, testFile)
    for i=1:length(files)
        testFile = files{i};
        currDir   = fileparts(which('testLegacy'));
        legacyDir = fullfile(currDir, 'legacy');

        % Clean up prior outputs
        seqFile = fullfile(legacyDir, [testFile '.seq']);
        outFile = fullfile(legacyDir, [testFile '.matlab.out']);
        if exist(seqFile, 'file'), delete(seqFile); end
        if exist(outFile, 'file'), delete(outFile); end

        % Run the legacy script from inside legacyDir
        od = cd(legacyDir);
        %cleanupObj = onCleanup(@() cd(od));
        try
            evalin('caller',testFile);
        catch ME
            testCase.verifyFail(...
                sprintf('Script %s threw: %s', testFile, ME.message));
            cd(od);
            return;
        end

        cd(od);

        % Compare output against approved reference
        approvedFile = fullfile(legacyDir, 'approved', [testFile '.matlab.out']);
        if mr.aux.isOctave()
            if exist(outFile, 'file') ~= 2
                error('Script %s did not produce %s.out', testFile, testFile);
            end
            cmpRes = compareTextFiles(outFile, approvedFile);
            if ~cmpRes
                error('Output of %s differs from approved reference', testFile);
            end
            fprintf('Succeded legacy integration test for sequence: %s\n', testFile);
        else
            testCase.verifyTrue(exist(outFile, 'file') == 2, ...
                sprintf('Script %s did not produce %s.out', testFile, testFile));
            cmpRes = compareTextFiles(outFile, approvedFile);
            testCase.verifyTrue(cmpRes, ...
                sprintf('Output of %s differs from approved reference', testFile));
            testCase.log(1, ['Succeded legacy integration test for sequence: ' testFile]);
        end
    end
end

