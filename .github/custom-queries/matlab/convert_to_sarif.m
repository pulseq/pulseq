% Load the results from the .mat file
load('code-analysis-results.mat', 'results');

% Initialize SARIF structure
sarif.version = "2.1.0";
sarif.runs = struct([]);

% Initialize run object
run.tool.driver.name = "Octave Static Code Analysis for MATLAB";
run.tool.driver.version = "1.0";
run.results = [];

% Populate SARIF results
for i = 1:length(results)
    filePath = results{i}{1};
    issues = results{i}{2};
    
    % Loop over each issue found in the file
    for j = 1:length(issues)
        % Create SARIF result entry for each issue
        sarifResult = struct();
        sarifResult.ruleId = "MATLABIssue"; % Identifier for the issue type
        sarifResult.level = "warning";
        sarifResult.message = struct('text', issues{j});
        
        % Specify the location
        sarifResult.locations = struct([]);
        sarifResult.locations(1).physicalLocation.artifactLocation.uri = filePath;
        sarifResult.locations(1).physicalLocation.region.startLine = 1; % Set startLine if available
        
        % Add result to the SARIF run
        run.results = [run.results, sarifResult];
    end
end

% Add the populated run to SARIF
sarif.runs = [sarif.runs, run];

% Save SARIF output to a JSON file
jsonText = jsonencode(sarif);
fid = fopen('code-analysis-results.sarif', 'w');
if fid == -1
    error('Could not create SARIF file');
end
fwrite(fid, jsonText, 'char');
fclose(fid);
