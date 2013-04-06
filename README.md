# conjugate-gradient
Solves sparse symmetric positive definite linear systems.  These problems arise in many physical applications, like [linear elasticity](http://en.wikipedia.org/wiki/Linear_elasticity), [heat transfer](http://en.wikipedia.org/wiki/Heat_transfer) and other diffusion based transport phenomena.

This code implements the [conjugate gradient method](http://en.wikipedia.org/wiki/Conjugate_gradient_method) using a [Jacobi preconditioner](http://en.wikipedia.org/wiki/Preconditioner).

## Install

    npm install conjugate-gradient

## Example

```javascript
var pcg = require("conjugate-gradient")
  , CSRMatrix = require("csr-matrix")

//Create a matrix
var A = CSRMatrix.fromDense([[-2, 1, 0],
                             [ 1,-2, 1],
                             [ 0, 1,-2]])

//Create input vector
var B = new Float64Array([1, 0, 0])

//Solve equation:
//
//  A x = B
//
console.log(pcg(A, b))
```

### `require("conjugate-gradient")(A, b[, x0, tolerance, max_iter])`
Solves the equation Ax = b by conjugate gradient

* `A` is a symmetric positive definite matrix represented as a [CSRMatrix](https://github.com/mikolalysenko/csr-matrix)
* `b` is an array of length n
* `x0` is an optional initial guess for the solution to the equation.  If specified, the result of the solution will also get stored in this array
* `tolerance` is a cutoff tolerance for the solution.  (Default is 1e-5)
* `max_iter` is the maximum number of iterations to run the solver.  (Default is min(n, 20))

**Returns** An array encoding the solution to the equation Ax = b

# Credits
(c) 2013 Mikola Lysenko. MIT License

