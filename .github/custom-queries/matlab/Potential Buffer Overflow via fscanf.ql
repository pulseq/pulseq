/**
 * @name Potential buffer overflow in fscanf
 * @description Detects potential buffer overflow vulnerabilities when using `fscanf` without checking input size.
 * @kind problem
 * @problem.severity critical
 * @tags security, buffer-overflow
 */
import javascript

class BufferOverflowFscanf extends Expr {
  BufferOverflowFscanf() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr fscanfCall |
      fscanfCall.getSource() = this and
      fscanfCall.toString().matches("fscanf%") and
      not exists(Expr sizeCheck | 
        sizeCheck.getSource() = fscanfCall and
        sizeCheck.toString().matches("(size|length)%"))
    )
  }
}

from BufferOverflowFscanf fscanfCall
select fscanfCall, "Potential buffer overflow detected in fscanf usage without input size validation."
