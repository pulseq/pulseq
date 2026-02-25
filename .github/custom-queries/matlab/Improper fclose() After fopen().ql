/**
 * @name Missing fclose after fopen
 * @description Detects cases where a file is opened but not closed, leading to resource leaks.
 * @kind problem
 * @problem.severity warning
 * @tags security, resource-leak
 */
import javascript

class MatlabMissingFclose extends Expr {
  MatlabMissingFclose() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr fopenCall |
      fopenCall.getSource() = this and
      fopenCall.toString().
