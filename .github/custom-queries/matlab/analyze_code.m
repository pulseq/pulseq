% Add the current working directory and subdirectories to the path
addpath(genpath(pwd));

% Find all .m files in the current directory and subdirectories
files = dir('**/*.m');

% Initialize results
results = {};

% Loop through each file and perform analysis
for k = 1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);
    fid = fopen(filePath, 'r');

    if fid == -1
        disp(['Could not open file: ', filePath]);
        continue;
    end

    code = fread(fid, '*char')'; % Read file contents
    fclose(fid);

    % Initialize issues for the current file
    issues = {};

    % Perform various analysis checks (as before)...

    % Example check: Check for improper input handling in mathematical operations
    if contains(code, '+') || contains(code, '-') || contains(code, '*') || contains(code, '/') || ...
        contains(code, '^') || contains(code, 'sqrt') || contains(code, 'log')
        if ~contains(code, 'validate') && ~contains(code, 'sanitize')
            issues{end+1} = struct('ruleId', 'ImproperInputHandling', 'message', 'Improper input handling in mathematical operations detected.');
        end
    end

    % Store file path and issues found in results
    if ~isempty(issues)
        results{end+1} = {filePath, issues};
    end
end

% Save results as a .mat file
save('code-analysis-results.mat', 'results');

% Save results as JSON file
fid = fopen('code-analysis-results.json', 'w');
fwrite(fid, jsonencode(results));
fclose(fid);
