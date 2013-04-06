var pcg = require("../pcg.js")
  , CSRMatrix = require("csr-matrix")
  , almostEqual = require("almost-equal")

require("tap").test("pcg", function(t) {

  function testMatrix(A, b, tol) {
    console.log(b)
    var x = pcg(A, b)
      , Ax = A.apply(x)
    for(var i=0; i<Ax.length; ++i) {
      t.assert(almostEqual(Ax[i], b[i], tol, tol), Ax[i] + " = " + b[i])
    }
  }

  //Create a matrix
  var A = CSRMatrix.fromDense([[-2, 1, 0],
                               [ 1,-2, 1],
                               [ 0, 1,-2]])

  //Create input vector
  var b = new Float64Array([1, 0, 0])

  testMatrix(A, b, 1e-5)


  t.end()
})