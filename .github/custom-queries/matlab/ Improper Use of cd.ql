/**
 * @name Improper use of cd command
 * @description Detects usage of `cd` to change directories without validating whether the directory exists, potentially leading to errors.
 * @kind problem
 * @problem.severity warning
 * @tags security, directory-management
 */
import javascript

class ImproperCdUsage extends Expr {
  ImproperCdUsage() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr cdCall |
      cdCall.getSource() = this and
      cdCall.toString().matches("cd%") and
      not exists(Expr check |
        check.getSource() = cdCall and check.toString().matches("(exist|isdir)%"))
    )
  }
}

from ImproperCdUsage cdCall
select cdCall, "Improper use of cd command without validating the directory path."
