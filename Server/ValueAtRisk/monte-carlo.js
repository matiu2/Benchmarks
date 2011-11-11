fs = require('fs');
MathExt = require('./MathExt.js');
FABRIC = require('Fabric').createClient();

var useFabric = false;
var prng = new MathExt.random.ExpGenerator;
prng.seed(0);

var numStocks = 10;
var numTradingDays = 252;
var dt = 1.0/numTradingDays;
var sqrtDT = Math.sqrt(dt);

var priceMeans = [];
for (var i=0; i<numStocks; ++i)
  priceMeans[i] = 25.0/numTradingDays;

var priceDevs = [];
for (var i=0; i<numStocks; ++i)
  priceDevs[i] = 25.0/numTradingDays;

var priceCorrelations = MathExt.randomCorrelation(numStocks, prng);
console.log("priceCorrelations:");
console.log(priceCorrelations);

var priceCovariance = [];
for (var i=0; i<numStocks; ++i) {
  priceCovariance[i] = [];
  for (var j=0; j<numStocks; ++j) {
    priceCovariance[i][j] = priceDevs[i] * priceDevs[j] * priceCorrelations[i][j];
  }
}
console.log("priceCovariance:");
console.log(priceCovariance);

var choleskyTrans = MathExt.choleskyTrans(priceCovariance);

var drifts = [];
for (var i=0; i<numStocks; ++i)
  drifts[i] = priceMeans[i] - priceCovariance[i][i]/2;

//var numTrials = 1;
//var numTrials = 256;
//var numTrials = 65536;
//var numTrials = 262144;
var numTrials = 1048576;

var valueAtRisk;
if (useFabric) {
  var params = FABRIC.DG.createNode("params");
  params.addMember('numTradingDays', 'Size', numTradingDays);
  params.addMember('dt', 'Float64', dt);
  params.addMember('sqrtDT', 'Float64', sqrtDT);
  params.addMember('choleskyTrans', 'Float64['+numStocks+']['+numStocks+']');
  params.setData('choleskyTrans', choleskyTrans);
  params.addMember('drifts', 'Float64['+numStocks+']');
  params.setData('drifts', drifts);

  var runTrialOp = FABRIC.DG.createOperator("runTrial");
  runTrial = fs.readFileSync('runTrial.kl', 'utf8').split('%NS%').join(numStocks);
  //console.log(runTrial);
  runTrialOp.setSourceCode('runTrial.kl', runTrial);
  runTrialOp.setEntryFunctionName('runTrial');
  if (runTrialOp.getDiagnostics().length > 0 ) {
    console.log(runTrialOp.getDiagnostics());
    throw "Compile errors, aborting";
  }

  var runTrialBinding = FABRIC.DG.createBinding();
  runTrialBinding.setOperator(runTrialOp);
  runTrialBinding.setParameterLayout([
    'self.index',
    'params.numTradingDays',
    'params.dt',
    'params.sqrtDT',
    'params.choleskyTrans',
    'params.drifts',
    'self.value'
  ]);

  var sortOp = FABRIC.DG.createOperator("sort");
  sort = fs.readFileSync('sort.kl', 'utf8').split('%NS%').join(numStocks);
  //console.log(sort);
  sortOp.setSourceCode('sort.kl', sort);
  sortOp.setEntryFunctionName('sort');
  if (sortOp.getDiagnostics().length > 0 ) {
    console.log(sortOp.getDiagnostics());
    throw "Compile errors, aborting";
  }

  var sortBinding = FABRIC.DG.createBinding();
  sortBinding.setOperator(sortOp);
  sortBinding.setParameterLayout([
    'self.value<>'
  ]);

  var trials = FABRIC.DG.createNode('trials');
  trials.setCount(numTrials);
  trials.setDependency(params, 'params');
  trials.addMember('value', 'Float64');
  trials.bindings.append(runTrialBinding);
  trials.bindings.append(sortBinding);
  if (trials.getErrors().length > 0) {
    console.log(trials.getErrors());
    throw "DG errors, aborting";
  }
  trials.evaluate();

  valueAtRisk = (numStocks * 100.0) - trials.getData('value', Math.round(numTrials*0.05));
}
else {
  trialResults = [];
  for (var trial=0; trial<numTrials; ++trial) {
    prng.seed(4096*(1+trial));

    //console.log("trial="+trial+" numTradingDays="+numTradingDays+" dt="+dt+" sqrtDT="+sqrtDT);
    //console.log("choleskyTrans="+choleskyTrans);
    //console.log("drifts="+drifts);
    var amounts = [];
    for (var i=0; i<numStocks; ++i)
      amounts[i] = 100;

    for (var day=1; day<=numTradingDays; ++day) {
      var Z = MathExt.random.normalVec(numStocks, prng);
      //console.log("Z = "+Z);
      var X = MathExt.mat.mulVec(choleskyTrans, Z);
      //console.log("X = "+X);
      for (var i=0; i<numStocks; ++i) {
        amounts[i] *= Math.exp(drifts[i]*dt + X[i]*sqrtDT);
      }
    }

    var value = 0.0;
    for (var i=0; i<numStocks; ++i)
      value += amounts[i];
    trialResults.push(value);
  }

  var sort = function (v) {
    var partition = function (a, begin, end, pivot) {
      var piv = a[pivot];
      a[pivot] = a[end-1];
      a[end-1] = piv;
      var store = begin;
      for (var i=begin; i<end-1; ++i) {
        if (a[i] <= piv) {
          var t = a[store];
          a[store] = a[i];
          a[i] = t;
          ++store;
        }
      }
      var t = a[end-1];
      a[end-1] = a[store];
      a[store] = t;
      return store;
    };

    var qsort = function (a, begin, end) {
      if (end - begin <= 1)
        return;
      else {
        var pivot = partition(a, begin, end, begin+Math.round((end-begin)/2));
        qsort(a, begin, pivot);
        qsort(a, pivot+1, end);
      }
    };

    return qsort(v, 0, v.length);
  };

  sort(trialResults);
  valueAtRisk = (numStocks * 100.0) - trialResults[Math.round(numTrials*0.05)];
}
console.log("ValueAtRisk = " + valueAtRisk);
