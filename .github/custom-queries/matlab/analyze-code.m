% Add the current working directory and subdirectories to the path
addpath(genpath(pwd));

% Find all .m files in the current directory and subdirectories
files = dir('**/*.m');

% Initialize SARIF results structure
sarifVersion = '2.1.0';
toolName = 'Octave Static Code Analyzer';
results = struct('version', sarifVersion, 'runs', [], 'tool', struct('driver', struct('name', toolName)));

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

    % Check for improper input handling in mathematical operations
    if contains(code, '+') || contains(code, '-') || contains(code, '*') || contains(code, '/') || ...
       contains(code, '^') || contains(code, 'sqrt') || contains(code, 'log')
        if ~contains(code, 'validate') && ~contains(code, 'sanitize')
            issues{end+1} = 'Improper input handling in mathematical operations detected.';
        end
    end

    % Check for improper use of cd command
    if contains(code, 'cd')
        if ~contains(code, 'exist') && ~contains(code, 'isdir')
            issues{end+1} = 'Improper use of cd command without validating the directory path.';
        end
    end

    % Check for unsafe mkdir usage
    if contains(code, 'mkdir')
        if ~contains(code, 'exist') && ~contains(code, 'isdir')
            issues{end+1} = 'Unsafe use of mkdir without checking for directory existence or permissions.';
        end
    end

    % Check for deprecated functions usage
    if contains(code, 'str2num') || contains(code, 'input') || contains(code, 'addpath')
        issues{end+1} = 'Deprecated function usage detected. Consider updating to supported alternatives.';
    end

    % Additional Checks
    % Check for hard-coded credentials, IPs, etc. (similar checks as before)

    % Store results in SARIF format if issues are found
    if ~isempty(issues)
        result = struct();
        result.locations = struct('physicalLocation', struct('artifactLocation', struct('uri', filePath)));
        result.annotations = struct('message', struct('text', issues), 'severity', 'warning');
        results.runs(end+1) = struct('tool', struct('driver', struct('name', toolName)), 'results', result);
    end
end

% Save results to a SARIF file
sarifFileName = 'code-analysis-results.sarif';
fid = fopen(sarifFileName, 'w');
fprintf(fid, '%s\n', jsonencode(results));
fclose(fid);
