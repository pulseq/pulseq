/**
 * @name Deprecated function usage
 * @description Detects the use of deprecated functions in MATLAB, which might introduce vulnerabilities or break in future releases.
 * @kind problem
 * @problem.severity warning
 * @tags security, deprecated
 */
import javascript

class DeprecatedFunctionUsage extends Expr {
  DeprecatedFunctionUsage() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr deprecatedCall |
      deprecatedCall.getSource() = this and
      deprecatedCall.toString().matches("(str2num|input|addpath)%")
    )
  }
}

from DeprecatedFunctionUsage deprecatedCall
select deprecatedCall, "Deprecated function usage detected. Consider updating to supported alternatives."
