/**
 * @name Insecure URL access without TLS
 * @description Detects `urlread` or `webread` function usage for URLs that don't enforce HTTPS, which can lead to insecure communication.
 * @kind problem
 * @problem.severity high
 * @tags security, communication
 */
import javascript

class InsecureURLAccess extends Expr {
  InsecureURLAccess() {
    this.getFile().getName().endsWith(".m") and
    exists(Expr urlReadCall |
      urlReadCall.getSource() = this and
      urlReadCall.toString().matches("(urlread|webread)%") and
      not urlReadCall.toString().matches("https%")
    )
  }
}

from InsecureURLAccess urlReadCall
select urlReadCall, "Insecure URL access detected, consider enforcing HTTPS."
