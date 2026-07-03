/**
 * @name Insecure use of load function
 * @description Detects use of the `load` function without specifying the file format, which can lead to security risks.
 * @kind problem
 * @problem.severity medium
 * @tags security, data-integrity
 */
import javascript

class InsecureLoadUsage extends Expr {
  InsecureLoadUsage() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr loadCall |
      loadCall.getSource() = this and
      loadCall.toString().matches("load%") and
      not loadCall.toString().matches("-ascii|-mat")
    )
  }
}

from InsecureLoadUsage loadCall
select loadCall, "Insecure use of the load function detected. Specify the file format to ensure data integrity."
