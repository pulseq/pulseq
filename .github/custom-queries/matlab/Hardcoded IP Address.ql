/**
 * @name Hardcoded IP address
 * @description Detects hardcoded IP addresses in the MATLAB code, which could indicate potential misconfigurations or security vulnerabilities.
 * @kind problem
 * @problem.severity medium
 * @tags security, misconfiguration
 */
import javascript

class HardcodedIP extends Expr {
  HardcodedIP() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr ipAddress |
      ipAddress.getSource() = this and
      ipAddress.toString().matches("(\\d{1,3}\\.){3}\\d{1,3}")
    )
  }
}

from HardcodedIP ipAddress
select ipAddress, "Hardcoded IP address detected, consider using configuration files instead."
