/**
 * @name Unchecked system command execution
 * @description Detects `system()` function usage with untrusted inputs, which can lead to command injection vulnerabilities.
 * @kind problem
 * @problem.severity critical
 * @tags security, command-injection
 */
import javascript

class UncheckedSystemCommand extends Expr {
  UncheckedSystemCommand() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr systemCall |
      systemCall.getSource() = this and
      systemCall.toString().matches("system%") and
      not exists(Expr validation |
        validation.getSource() = systemCall and
        validation.toString().matches("(validate|sanitize)%"))
    )
  }
}

from UncheckedSystemCommand systemCall
select systemCall, "Unchecked system command execution detected. Ensure input is properly validated or sanitized."
