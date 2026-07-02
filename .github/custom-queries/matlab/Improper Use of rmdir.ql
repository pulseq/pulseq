/**
 * @name Improper rmdir usage
 * @description Detects `rmdir` calls without error handling or validation.
 * @kind problem
 * @problem.severity high
 * @tags security, resource-management
 */
import javascript

class ImproperRmdir extends Expr {
  ImproperRmdir() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr rmdirCall |
      rmdirCall.getSource() = this and
      rmdirCall.toString().matches("rmdir%") and
      not exists(Expr errorCheck | 
        errorCheck.getSource() = rmdirCall and errorCheck.toString().matches("(exist|isdir)%"))
    )
  }
}

from ImproperRmdir rmdirCall
select rmdirCall, "Improper use of rmdir without error handling or path validation."
