classdef testLegacy < matlab.unittest.TestCase
% testLegacy  Parameterised regression tests for legacy scripts.
%
%   Enumerates *.out files in tests/legacy/approved/ and creates one test
%   point per file.  Each test runs the corresponding script, then compares
%   the freshly generated .out file against the approved reference, ignoring
%   line-ending differences.
%
%   Run with:   runtests('testLegacy')

    properties (TestParameter)
        testFile = testLegacy.enumerateFiles();
    end

    methods (Static)
        function files = enumerateFiles()
            currDir  = fileparts(which('testLegacy'));
            % test whether the main toolbox is in the path
            try
                sys=mr.opts();
            catch
                addpath(fullfile(currDir,'../matlab'));
            end
            approvedDir = fullfile(currDir, 'legacy', 'approved');
            d = dir(fullfile(approvedDir, '*.matlab.out'));
            files = struct();
            for k = 1:numel(d)
                [~, name, ~] = fileparts(d(k).name);
                if strcmp(name(end-6:end),'.matlab'), name=name(1:end-7); end
                safeName = matlab.lang.makeValidName(name);
                files.(safeName) = name;
            end
        end
    end

    methods (Test)
        function verifyLegacyFile(testCase, testFile)
            currDir   = fileparts(which('testLegacy'));
            legacyDir = fullfile(currDir, 'legacy');

            % Clean up prior outputs
            seqFile = fullfile(legacyDir, [testFile '.seq']);
            outFile = fullfile(legacyDir, [testFile '.matlab.out']);
            if exist(seqFile, 'file'), delete(seqFile); end
            if exist(outFile, 'file'), delete(outFile); end

            % Run the legacy script from inside legacyDir
            od = cd(legacyDir);
            cleanupObj = onCleanup(@() cd(od));
            try
                eval(testFile);
            catch ME
                testCase.verifyFail(...
                    sprintf('Script %s threw: %s', testFile, ME.message));
                return;
            end

            % Compare output against approved reference
            approvedFile = fullfile(legacyDir, 'approved', [testFile '.matlab.out']);
            testCase.verifyTrue(exist(outFile, 'file') == 2, ...
                sprintf('Script %s did not produce %s.out', testFile, testFile));
            cmpRes = compareTextFiles(outFile, approvedFile);
            testCase.verifyTrue(cmpRes, ...
                sprintf('Output of %s differs from approved reference', testFile));
        end
    end
end

% ---- local helper (accessible within classdef file) --------------------
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
