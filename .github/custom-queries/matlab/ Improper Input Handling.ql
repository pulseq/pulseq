/**
 * @name Improper input handling in mathematical operations
 * @description Detects mathematical operations on untrusted inputs without proper validation or sanitization, leading to potential data integrity issues.
 * @kind problem
 * @problem.severity high
 * @tags security, input-validation
 */
import javascript

class ImproperInputInMathOp extends Expr {
  ImproperInputInMathOp() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr mathOp |
      mathOp.getSource() = this and
      mathOp.toString().matches("(\\+|\\-|\\*|\\/|\\^|sqrt|log)%") and
      not exists(Expr inputValidation |
        inputValidation.getSource() = mathOp and inputValidation.toString().matches("(validate|sanitize)%"))
    )
  }
}

from ImproperInputInMathOp mathOp
select mathOp, "Improper input handling in mathematical operations detected. Ensure inputs are validated or sanitized."
