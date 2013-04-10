"use strict"

var bits = require("bits")
  , almostEqual = require("almost-equal")

var R = new Float64Array(1024)
var P = new Float64Array(1024)
var D = new Float64Array(1024)
var Z = new Float64Array(1024)

function reserve(n) {
  if(n < R.length) {
    return
  }
  var nsize = bits.nextPow2(n)
  R = new Float64Array(nsize)
  P = new Float64Array(nsize)
  Z = new Float64Array(nsize)
  D = new Float64Array(nsize)
}

function conjugateGradient(A, b, x, tolerance, max_iter) {
  var abs = Math.abs
    , max = Math.max
    , EPSILON = almostEqual.FLT_EPSILON
    , n = A.rowCount
    , i, j, k
    , alpha_n, alpha_d, alpha, beta, rnorm, s
  if(!tolerance) {
    tolerance = 1e-5
  }
  if(!max_iter) {
    max_iter = Math.min(n, 20)
  }
  if(!x) {
    if(b.buffer) {
      x = new b.constructor(b.buffer.slice(0))
    } else {
      x = b.slice(0)
    }
  }
  reserve(n)
  //Compute preconditioner
  for(i=0; i<n; ++i) {
    s = A.get(i, i)
    if(abs(s) > EPSILON) {
      D[i] = 1.0 / s
    } else {
      D[i] = 1.0
    }
  }
  //Initialize 
  A.apply(x, R)
  for(i=0; i<n; ++i) {
    R[i] = b[i] - R[i]
    Z[i] = D[i] * R[i]
    P[i] = Z[i]
  }
  //Iterate
  for(k=0; k<max_iter; ++k) {
    alpha_n = 0.0
    for(i=0; i<n; ++i) {
      alpha_n += R[i] * Z[i]
    }
    A.apply(P, Z)
    alpha_d = 0.0
    for(i=0; i<n; ++i) {
      alpha_d +=  P[i] * Z[i]
    }
    alpha = alpha_n / alpha_d
    beta = 0.0
    rnorm = 0.0
    for(i=0; i<n; ++i) {
      x[i] += alpha * P[i]
      R[i] -= alpha * Z[i]
      Z[i]  = D[i] * R[i]
      beta += R[i] * Z[i]
      rnorm = max(rnorm, abs(R[i]))
    }
    if(rnorm < tolerance) {
      break
    }
    beta /= alpha_n
    for(i=0; i<n; ++i) {
      P[i] = Z[i] + beta * P[i]
    }
  }
  return x
}

module.exports = conjugateGradient