/**
 * @name Untrusted input in eval
 * @description Detects `eval` function usage with untrusted or unsanitized input, leading to potential code injection.
 * @kind problem
 * @problem.severity critical
 * @tags security, code-injection
 */
import javascript

class EvalUntrustedInput extends Expr {
  EvalUntrustedInput() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr evalCall |
      evalCall.getSource() = this and
      evalCall.toString().matches("eval%") and
      not exists(Expr inputValidation |
        inputValidation.getSource() = evalCall and inputValidation.toString().matches("(validate|sanitiz)%"))
    )
  }
}

from EvalUntrustedInput evalCall
select evalCall, "Potential code injection via untrusted input passed to eval function."
