/**
 * @name Hardcoded credentials
 * @description Detects potential hardcoded credentials, such as passwords or API keys, in MATLAB scripts.
 * @kind problem
 * @problem.severity critical
 * @tags security, hardcoded-credentials
 */

import javascript

class MatlabHardcodedCredentials extends Expr {
  MatlabHardcodedCredentials() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr hardcoded |
      hardcoded.getSource() = this and
      hardcoded.toString().matches("(password|passwd|apiKey|secret|token)\\s*=\\s*['\"]\\w+['\"]")
    )
  }
}

from MatlabHardcodedCredentials hardcoded
select hardcoded, "Potential hardcoded credentials detected in MATLAB code."
