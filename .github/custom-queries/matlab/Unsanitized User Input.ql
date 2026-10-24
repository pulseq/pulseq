/**
 * @name Unsanitized user input in plotting functions
 * @description Detects unsanitized or unvalidated user input passed to plotting functions, which can lead to vulnerabilities.
 * @kind problem
 * @problem.severity medium
 * @tags security, input-validation
 */
import javascript

class UnsanitizedInputInPlot extends Expr {
  UnsanitizedInputInPlot() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr plotCall |
      plotCall.getSource() = this and
      plotCall.toString().matches("(plot|scatter|bar|hist)%") and
      not exists(Expr inputValidation |
        inputValidation.getSource() = plotCall and
        inputValidation.toString().matches("(validate|sanitize)%"))
    )
  }
}

from UnsanitizedInputInPlot plotCall
select plotCall, "Unsanitized user input detected in plotting function. Ensure input is properly validated or sanitized."
