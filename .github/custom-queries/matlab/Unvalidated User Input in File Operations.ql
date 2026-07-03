/**
 * @name Unvalidated input in file operations
 * @description Detects unvalidated user input used in file operations, leading to potential file injection or manipulation.
 * @kind problem
 * @problem.severity critical
 * @tags security, input-validation
 */
import javascript

class UnvalidatedInputFileOp extends Expr {
  UnvalidatedInputFileOp() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr fileOp |
      fileOp.getSource() = this and
      fileOp.toString().matches("(fopen|fwrite|fclose)%") and
      not exists(Expr inputValidation |
        inputValidation.getSource() = fileOp and inputValidation.toString().matches("(validate|sanitiz)%"))
    )
  }
}

from UnvalidatedInputFileOp fileOp
select fileOp, "Unvalidated input passed to file operations like fopen or fwrite, leading to security risks."
