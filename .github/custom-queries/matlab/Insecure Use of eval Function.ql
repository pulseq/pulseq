/**
 * @name Insecure use of eval
 * @description Detects insecure usage of the `eval` function, which can lead to code injection vulnerabilities.
 * @kind problem
 * @problem.severity critical
 * @tags security, code-injection
 */

import javascript

class MatlabInsecureEvalUsage extends Expr {
  MatlabInsecureEvalUsage() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr evalCall |
      evalCall.getSource() = this and
      evalCall.toString().matches("eval%")
    )
  }
}

from MatlabInsecureEvalUsage evalCall
select evalCall, "Insecure use of eval function, consider alternatives or validate input."
