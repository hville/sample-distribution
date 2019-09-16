var upperBound = require('./src/upperbound')

module.exports = CDF

function CDF(len) {
	// properties
	this.vs = Array(len)
	this.rs = []
}

CDF.prototype = {
	constructor: CDF,
	get N() { return this.rs[this.rs.length-1] },
	get E() {
		var vs = this.vs,
				rs = this.rs,
				M1 = rs.length-1,
				sum = vs[0] + vs[M1]
		for (var i=0; i<M1; ++i) sum += (vs[i+1] + vs[i]) * (rs[i+1] - rs[i])
		return sum/2/rs[M1]
	},
	R: function(r) { return this.Q( (r-0.5)/this.rs[this.rs.length-1] ) },
	Q: quantile,
	F: cdf,
	f: pdf,
	push: push
}

function push(x) {
	var j = upperBound(this.vs, x),
			vs = this.vs,
			rs = this.rs,
			M = rs.length,
			ir = 0
	// lossless insert
	if (M < vs.length) {
		for (ir=rs.length; ir>j; --ir) {
			rs[ir] = ir+1
			vs[ir] = vs[ir-1]
		}
		rs[j] = j ? rs[j-1] + 1 : 1
		vs[j] = x
	}
	// compression, droping a value while maintaining the interval average
	else if (j === M) newmax(vs, rs, x)
	else if (j === 0) {
		var u = vs[0],
				Δwx = vs[2] - x
		if (Δwx !== 0) {
			//KEEP 0 : r'(w-x) = (w+v-2x) + r(w-u) - s(v-u)
			//DROP 0 : r'(w-x) = (w+v-2x) + r(w-u)
			var r_v = (vs[2] + vs[1] - x - x + rs[1]*(vs[2] - u)) / Δwx,
					r_u = r_v - rs[2] * (vs[1]-u) / Δwx
			if ( (u + vs[1] > x + vs[2] && r_u > rs[0]) || r_v > rs[2]+1) {
				vs[1] = u
				rs[1] = r_u
			}
			else rs[1] = r_v
			vs[0] = x
		}
		for (ir=2; ir<rs.length; ++ir) ++rs[ir]
	}
	else {
		var i = j === 1 || (j !== M-1 && x+x > vs[j]+vs[j-1] ) ? j : j-1
		var w = vs[i+1],
				v = vs[i],
				Δwu = w-vs[i-1]
		if (Δwu !== 0) {
			var r_x = rs[i] + ( w-x + (x-v)*(rs[i+1] - rs[i-1]) ) / Δwu
			if ( v < x
				? vs[i-1]+w < v+x || r_x > rs[i+1]+1
				: v+x < vs[i-1] + w || r_x < rs[i-1]
			) rs[i] += ( w-x + v-x ) / Δwu
			else {
				rs[i] = r_x
				vs[i] = x
			}
		}
		while(++i<rs.length) ++rs[i]
	}
}

/**
 * Quantile function, provide the value for a given probability
 * @param {number} prob - probability or array of probabilities
 * @return {number} value or array of values
 */
function quantile(prob) {
	var vs = this.vs,
			rs = this.rs,
			M = rs.length,
			h = rs[M-1] * prob + 0.5, // 0.5 <= h <= N+0.5
			j = upperBound(rs, h), //      0 <= j <= M
			i = j-1
	return j === 0 ? vs[0]
		: j === M ? vs[M-1]
		: vs[i] + (vs[j] - vs[i]) * (h-rs[i]) / (rs[j]-rs[i])
}


/**
 * @param {number} x - probability or array of probabilities
 * @return {number} value or array of values
 */
function cdf(x) {
	var vs = this.vs,
			rs = this.rs,
			M = rs.length,
			N = rs[M-1],
			j = upperBound(vs, x),
			i = j-1
	return (j === 0 ? 0.5
		: j === M ? (N - 0.5)
		: rs[i] - 0.5 + (rs[j] - rs[i]) * (x - vs[i]) / (vs[j] - vs[i])
	) / N
}

/**
 * @param {number} x - probability or array of probabilities
 * @return {number} value or array of values
 */
function pdf(x) {
	var vs = this.vs,
			rs = this.rs,
			M = rs.length,
			N = rs[M-1]
	if (x === vs[0] || x === vs[M-1]) return 0.5/N
	var j = upperBound(vs, x),
			i = j-1
	return j === 0 || j === M ? 0 : (rs[j] - rs[i]) / N
}

function newmax(vs, rs, x) {
	var j = rs.length-1,
			i = j-1,
			h = i-1,
			Δwv = vs[j]-vs[i],
			Δxu = x - vs[h],
			rjh = rs[i]*(vs[j]-vs[h])
	if (Δxu !== 0) {
		var r_w = ( rs[j]*(x-vs[i]) + rjh - rs[h]*Δwv ) / Δxu,
				r_v = ( rs[j]*(x-vs[j]) + rjh - Δwv ) / Δxu
		if ( r_v < rs[h] || (vs[i] + vs[j] < vs[h] + x && r_w < rs[j]+1)) {
			vs[i] = vs[j]
			rs[i] = r_w
		}
		else rs[i] = r_v
		vs[j] = x
	}
	++rs[j]
}
