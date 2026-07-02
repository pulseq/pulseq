/**
 * @name Unsafe fopen usage
 * @description Detects `fopen` calls without error handling or proper file permissions.
 * @kind problem
 * @problem.severity warning
 * @tags security, file-handling
 */

import javascript

class MatlabUnsafeFopenUsage extends Expr {
  MatlabUnsafeFopenUsage() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr fopenCall |
      fopenCall.getSource() = this and
      fopenCall.toString().matches("fopen%") and
      not exists(Expr fcloseCall | 
        fcloseCall.getSource() = fopenCall and 
        fcloseCall.toString().matches("fclose%"))
    )
  }
}

from MatlabUnsafeFopenUsage fopenCall
select fopenCall, "File opened with fopen but no fclose or error handling detected."
