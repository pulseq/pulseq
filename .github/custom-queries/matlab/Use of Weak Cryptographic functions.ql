/**
 * @name Weak cryptographic functions
 * @description Detects usage of weak or deprecated cryptographic functions such as MD5.
 * @kind problem
 * @problem.severity high
 * @tags security, cryptography
 */
import javascript

class WeakCryptography extends Expr {
  WeakCryptography() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr cryptoCall |
      cryptoCall.getSource() = this and
      cryptoCall.toString().matches("(md5|sha1)%")
    )
  }
}

from WeakCryptography cryptoCall
select cryptoCall, "Weak cryptographic function detected (md5 or sha1), consider stronger alternatives."
