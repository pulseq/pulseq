/**
 * @name Insecure file permissions
 * @description Detects insecure file permission changes through `chmod` or `chown` commands.
 * @kind problem
 * @problem.severity high
 * @tags security, file-permissions
 */

import javascript

class MatlabInsecureFilePermissions extends Expr {
  MatlabInsecureFilePermissions() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr permCall |
      permCall.getSource() = this and
      permCall.toString().matches("(chmod|chown)%")
    )
  }
}

from MatlabInsecureFilePermissions permCall
select permCall, "Insecure file permission changes detected, review chmod/chown usage."