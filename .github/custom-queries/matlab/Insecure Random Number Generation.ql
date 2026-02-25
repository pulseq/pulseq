/**
 * @name Insecure random number generation
 * @description Detects usage of non-cryptographically secure random number generators like `rand` and `randi`.
 * @kind problem
 * @problem.severity warning
 * @tags security, crypto
 */

import javascript

class MatlabInsecureRandUsage extends Expr {
  MatlabInsecureRandUsage() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr randCall |
      randCall.getSource() = this and
      randCall.toString().matches("rand%")
    )
  }
}

from MatlabInsecureRandUsage randCall
select randCall, "Insecure random number generator used (rand or randi), consider using a secure alternative."
