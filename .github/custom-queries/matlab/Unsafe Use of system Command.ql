/**
 * @name Insecure use of system command
 * @description Detects unvalidated input passed to the `system` function, which can lead to command injection.
 * @kind problem
 * @problem.severity critical
 * @tags security, command-injection
 */

import javascript

class MatlabSystemCommand extends Expr {
  MatlabSystemCommand() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr sysCall |
      sysCall.getSource() = this and
      sysCall.toString().matches("system%") and
      not exists(Expr inputValidation |
        inputValidation.getSource() = sysCall and
        inputValidation.toString().matches("validate%"))
    )
  }
}

from MatlabSystemCommand sysCall
select sysCall, "Potentially unsafe system command execution without input validation."
