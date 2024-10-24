% Add the current working directory and subdirectories to the path
addpath(genpath(pwd));

% Find all .m files in the current directory and subdirectories
files = dir('**/*.m');

% Initialize results
results = struct('version', '2.1.0', 'runs', []);

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

    % 1. Check for hard-coded credentials
    if contains(code, 'password') || contains(code, 'passwd') || contains(code, 'apiKey') || ...
       contains(code, 'secret') || contains(code, 'token')
        if ~contains(code, 'validate') && ~contains(code, 'sanitize')
            issues{end+1} = 'Hard-coded credentials detected. Avoid hard-coding sensitive information.';
        end
    end

    % 2. Check for hard-coded IP addresses
    if ~isempty(regexp(code, '(\d{1,3}\.){3}\d{1,3}', 'once')) % Regular expression for matching IP addresses
        issues{end+1} = 'Hard-coded IP address detected. Consider using a configuration file or environment variables.';
    end

    % 3. Check for improper use of pause function
    if contains(code, 'pause')
        if ~contains(code, 'check') && ~contains(code, 'validate')
            issues{end+1} = 'Improper use of pause function detected without timing checks.';
        end
    end

    % 4. Check for improper use of rmdir
    if contains(code, 'rmdir')
        if ~contains(code, 'exist') && ~contains(code, 'isdir')
            issues{end+1} = 'Improper use of rmdir without validating the directory path.';
        end
    end

    % 5. Check for missing fclose after fopen
    if contains(code, 'fopen')
        if ~contains(code, 'fclose')
            issues{end+1} = 'Missing fclose after fopen. Ensure to close files to prevent resource leaks.';
        end
    end

    % 6. Check for insecure file permissions
    if contains(code, 'chmod') || contains(code, 'chown')
        issues{end+1} = 'Insecure file permissions command detected. Review file permission handling.';
    end

    % 7. Check for insecure random number generation
    if contains(code, 'rand')
        issues{end+1} = 'Insecure random number generation detected. Consider using a cryptographically secure alternative.';
    end

    % 8. Check for insecure save function usage
    if contains(code, 'save') && ~contains(code, '-v7.3')
        issues{end+1} = 'Insecure save function usage detected. Consider using the -v7.3 option for saving files.';
    end

    % 9. Check for insecure use of eval
    if contains(code, 'eval')
        issues{end+1} = 'Potentially unsafe use of eval detected. Avoid using eval for executing arbitrary code.';
    end

    % 10. Check for unclosed figures
    if contains(code, 'figure') && ~contains(code, 'close')
        issues{end+1} = 'Figures created but not closed. Ensure to close figures to free up system resources.';
    end

    % 11. Check for lack of error handling
    if contains(code, 'try') && ~contains(code, 'catch')
        issues{end+1} = 'Try without catch. Ensure proper error handling is implemented.';
    end

    % 12. Check for usage of global variables
    if contains(code, 'global')
        issues{end+1} = 'Use of global variables detected. Minimize the use of global variables for better encapsulation.';
    end

    % 13. Check for comments not matching code
    if contains(code, '%')
        % Check for comment length versus code length
        lines = strsplit(code, '\n');
        for line = lines
            if length(line{1}) < 20 && contains(line{1}, '%') % Short comments
                issues{end+1} = 'Comments should provide sufficient context; consider expanding short comments.';
            end
        end
    end

    % 14. Check for excessively long lines
    lines = strsplit(code, '\n');
    for i = 1:length(lines)
        if length(lines{i}) > 80
            issues{end+1} = sprintf('Line %d exceeds 80 characters. Consider breaking it into multiple lines.', i);
        end
    end

    % Store file path and issues found in results
    if ~isempty(issues)
        issueInstances = struct('ruleId', 'CustomRule', 'message', issues, 'locations', struct('physicalLocation', struct('artifactLocation', struct('uri', filePath))));
        results.runs(end+1).results(end+1) = issueInstances;
    end
end

% Save results to SARIF format
sarifFileName = 'code-analysis-results.sarif';
fid = fopen(sarifFileName, 'w');
fprintf(fid, '{"version": "2.1.0", "runs": [%s]}', jsonencode(results.runs));
fclose(fid);
disp(['Results saved to ', sarifFileName]);
