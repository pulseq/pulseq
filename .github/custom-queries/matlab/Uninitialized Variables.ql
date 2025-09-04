/**
 * @name Uninitialized variable in MATLAB
 * @description Find variables used without being initialized.
 * @kind problem
 * @problem.severity warning
 * @tags reliability
 */

import javascript

// Find a variable usage before its assignment
class MatlabVariableUsage extends Expr {
  MatlabVariableUsage() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr usage |
      usage = this and
      not exists(Expr init | 
        init.getSource() = usage and
        init.getFile().getName().endsWith(".m"))
    )
  }
}

from MatlabVariableUsage usage
select usage, "Possible uninitialized variable usage in MATLAB"
