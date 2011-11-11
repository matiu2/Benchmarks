module.exports = {
  random: {
    JSBuiltinGenerator: function () {
      this.seed = function (index) {
      };
      this.nextUniform01 = function () {
        return Math.random();
      };
    },

    ExpGenerator: function () {
      var base = { lo: 300773&0x7FFF, hi: 300773>>15 };
      var mul = function (lhs, rhs) {
        return {
          lo: (lhs.lo * rhs.lo)&0x7FFF,
          hi: (
              ((lhs.lo*rhs.lo)>>15)
              + lhs.lo * rhs.hi
              + lhs.hi * rhs.lo
              )&0x7FFF
        };
      };
      var cur;
      this.seed = function (index) {
        cur = { lo: 1, hi: 0 };
        var mask = 1;
        var mult = { lo: base.lo, hi: base.hi };
        for (var i=0; i<30; ++i) {
          if ( index & mask )
            cur = mul(cur, mult);
          mult = mul(mult, mult);
          mask <<= 1;
        }
      };
      this.nextUniform01 = function () {
        cur = mul(cur, base);
        var seed = ((cur.lo + (cur.hi<<15))&0x3FFFFFFF);
        //console.log("seed = " + seed);
        return seed / 1073741824.0;
      }
    },

    normal: function (generator) {
      var x1, x2, w;
      do {
        x1 = 2.0 * generator.nextUniform01() - 1.0;
        x2 = 2.0 * generator.nextUniform01() - 1.0;
        w = x1 * x1 + x2 * x2;
      } while (w >= 1.0);
      w = Math.sqrt( (-2.0 * Math.log(w)) / w );
      return x1 * w;
    },

    normalVec: function (n, generator) {
      var v = [];
      for (var i=0; i<n; ++i)
        v[i] = this.normal(generator);
      return v;
    }
  },

  mat: {
    zero: function (rows, cols) {
      var M = [];
      for (var i=0; i<rows; ++i) {
        M[i] = [];
        for (var j=0; j<cols; ++j)
          M[i][j] = 0.0;
      }
      return M;
    },

    mulVec: function (M, v) {
      var r = [];
      for (var i=0; i<M.length; ++i) {
        var sum = 0.0;
        for (var j=0; j<v.length; ++j)
          sum += M[i][j] * v[j];
        r[i] = sum;
      }
      return r;
    },

    mul: function (M, N) {
      var R = [];
      for (var i=0; i<M.length; ++i) {
        if (M[i].length != N.length)
          throw "bad matricies";
        R[i] = [];
        for (var j=0; j<N[0].length; ++j) {
          var sum = 0;
          for (var k=0; k<N.length; ++k)
            sum += M[i][k] * N[k][j];
          R[i][j] = sum;
        }
      }
      return R;
    },

    trans: function (M) {
      var R = [];
      for (var i=0; i<M[0].length; ++i) {
        R[i] = [];
        for (var j=0; j<M.length; ++j) {
          R[i][j] = M[j][i];
        }
      }
      return R;
    }
  },

  randomCorrelation: function (n, prng) {
    var T = [];

    // Generate n uniform random columns
    for (var i=0; i<n; ++i) {
      T[i] = [];
      for (var j=0; j<n; ++j) {
        T[i][j] = this.random.normal(prng);
      }
    }

    // Normalize the columns
    for (var j=0; j<n; ++j) {
      var sqSum = 0;
      for (var i=0; i<n; ++i)
        sqSum += T[i][j] * T[i][j];
      var norm = Math.sqrt(sqSum);
      for (var i=0; i<n; ++i)
        T[i][j] /= norm;
    }

    // result is T'*T
    return this.mat.mul(this.mat.trans(T), T);
  },

  choleskyTrans: function (A) {
    var n = A.length;
    
    var L = this.mat.zero(n, n);
    for (var i=0; i<n; ++i) {
      for (var j=0; j<(i+1); ++j ) {
        var s = 0.0;
        for (var k=0; k<j; ++k)
          s += L[i][k] * L[j][k];
        if (i == j)
          L[i][i] = Math.sqrt(A[i][i]-s);
        else
          L[i][j] = 1.0/L[j][j]*(A[i][j]-s);
      }
    }
    return L;
  },

  cholesky: function (A) {
    return this.mat.trans(this.choleskyTrans(A));
  }
};
