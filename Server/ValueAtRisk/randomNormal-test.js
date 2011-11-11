var MathExt = require('./MathExt.js');

var bins = [];
for (var j=0; j<=100; ++j)
  bins[j] = 0;

for (var i=0; i<10000; ++i) {
  var randomNormal = MathExt.randomNormal();
  var adjRanNorm = Math.round(randomNormal*10 + 50);
  if (adjRanNorm >= 0 && adjRanNorm <= 100)
    ++bins[adjRanNorm];
}

var stars = "****************************************************************************************************************************************************************************************************************************";
for (var j=0; j<=100; ++j) {
  console.log(stars.substr(0, bins[j]));
}
