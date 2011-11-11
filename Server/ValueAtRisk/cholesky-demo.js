var MathExt = require('./MathExt.js');

var test = function(M) {
  console.log("M:");
  console.log(M);
  console.log("MathExt.cholesky(M)");
  console.log(MathExt.cholesky(M));
};

test([
  [25, 15, -5],
  [15, 18,  0],
  [-5,  0, 11]
]);

test([
  [18, 22, 54, 42],
  [22, 70, 86, 62],
  [54, 86, 174, 134],
  [42, 62, 134, 106]
]);
