/**
 * @name Unsafe mkdir usage
 * @description Detects the use of `mkdir` without checking for the directory's existence or permissions.
 * @kind problem
 * @problem.severity warning
 * @tags security, resource-management
 */
import javascript

class UnsafeMkdir extends Expr {
  UnsafeMkdir() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr mkdirCall |
      mkdirCall.getSource() = this and
      mkdirCall.toString().matches("mkdir%") and
      not exists(Expr check | 
        check.getSource() = mkdirCall and check.toString().matches("(exist|isdir)%"))
    )
  }
}

from UnsafeMkdir mkdirCall
select mkdirCall, "Unsafe use of mkdir without checking for directory existence or permissions."
