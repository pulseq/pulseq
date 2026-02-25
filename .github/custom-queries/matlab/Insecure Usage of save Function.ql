/**
 * @name Insecure save function usage
 * @description Detects insecure usage of the `save` function without specifying version or compression, leading to data risks.
 * @kind problem
 * @problem.severity warning
 * @tags security, data-integrity
 */
import javascript

class InsecureSaveFunction extends Expr {
  InsecureSaveFunction() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr saveCall |
      saveCall.getSource() = this and
      saveCall.toString().matches("save%") and
      not saveCall.toString().matches("-v7\\.3%")
    )
  }
}

from InsecureSaveFunction saveCall
select saveCall, "Insecure save function usage detected, consider using '-v7.3' for better data integrity."
