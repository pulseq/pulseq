/**
 * @name Improper use of pause function
 * @description Detects improper usage of the `pause` function, which can lead to inefficient code execution or resource misuse.
 * @kind problem
 * @problem.severity medium
 * @tags security, resource-management
 */
import javascript

class ImproperPauseUsage extends Expr {
  ImproperPauseUsage() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr pauseCall |
      pauseCall.getSource() = this and
      pauseCall.toString().matches("pause%") and
      not exists(Expr timingCheck |
        timingCheck.getSource() = pauseCall and timingCheck.toString().matches("(check|validate)%"))
    )
  }
}

from ImproperPauseUsage pauseCall
select pauseCall, "Improper use of the pause function detected, which could lead to resource misuse or inefficient code execution."
