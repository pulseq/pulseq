import scipy.io
import json

# Load the .mat file
data = scipy.io.loadmat('code-analysis-results.mat')
results = data.get('results', [])

# SARIF template
sarif = {
    "version": "2.1.0",
    "runs": [
        {
            "tool": {
                "driver": {
                    "name": "MATLAB Code Analysis",
                    "informationUri": "https://github.com/username/repo",
                    "rules": [
                        {
                            "id": "MATLAB001",
                            "name": "MATLAB Code Issue",
                            "shortDescription": {"text": "Identified MATLAB Code Issue"},
                            "fullDescription": {"text": "This issue was identified in the MATLAB code."},
                            "defaultConfiguration": {"level": "warning"}
                        }
                    ]
                }
            },
            "results": []
        }
    ]
}

# Convert each entry in the results to SARIF format
for item in results:
    file_path = item[0][0]
    issues = item[1][0]
    for issue in issues:
        sarif_issue = {
            "ruleId": "MATLAB001",
            "message": {"text": issue[0]},
            "locations": [
                {
                    "physicalLocation": {
                        "artifactLocation": {"uri": file_path},
                    }
                }
            ]
        }
        sarif["runs"][0]["results"].append(sarif_issue)

# Write to SARIF file
with open('code-analysis-results.sarif', 'w') as f:
    json.dump(sarif, f, indent=2)
