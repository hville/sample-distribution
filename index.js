module.exports = CDF

function CDF(len) {
	var buffer = len.constructor === len ? len : new ArrayBuffer(len << 7)
	this.vs = new Float64Array(buffer, 0, len)
	this.rs = new Float64Array(buffer, len << 6, len)
}

CDF.prototype = {
	constructor: CDF,
	// Number of samples
	get N() { return this.rs[this.rs.length-1] },
	// Expected Value
	get E() {
		var vs = this.vs,
				rs = this.rs,
				M = Math.min(rs.length, rs[rs.length-1]),
				Mm = M-1,
				sum = vs[0] + vs[Mm]
		for (var i=0; i<Mm; ++i) sum += (vs[i+1] + vs[i]) * (rs[i+1] - rs[i])
		return sum/2/rs[Mm]
	},
	// Origin Moments
	M: function(order) {
		var vs = this.vs,
				rs = this.rs,
				M = Math.min(rs.length, rs[rs.length-1]),
				Mm = M-1,
				Op = order + 1,
				sum = ( Math.pow(vs[0], order) + Math.pow(vs[Mm], order) ) * Op / 2
		for (var i=0; i<Mm; ++i) sum += (rs[i+1] - rs[i]) * ( Math.pow(vs[i+1], Op) - Math.pow(vs[i], Op) ) / (vs[i+1] - vs[i])
		return sum / Op / rs[Mm]
	},
	/**
	 * Quantile function, provide the value for a given probability
	 * @param {number} prob - probability or array of probabilities
	 * @return {number} value or array of values
	 */
	Q: function(prob) {
		var vs = this.vs,
				rs = this.rs,
				M = Math.min(rs.length, rs[rs.length-1]),
				h = rs[M-1] * prob + 0.5, // 0.5 <= h <= N+0.5
				j = topIndex(rs, h, M), //      0 <= j <= M
				i = j-1
		return j === 0 ? vs[0]
			: j === M ? vs[M-1]
			: vs[i] + (vs[j] - vs[i]) * (h-rs[i]) / (rs[j]-rs[i])
	},
	/**
	 * @param {number} x - probability or array of probabilities
	 * @return {number} value or array of values
	 */
	F: function(x) {
		var vs = this.vs,
				rs = this.rs,
				M = Math.min(rs.length, rs[rs.length-1]),
				N = rs[M-1],
				j = topIndex(vs, x, M),
				i = j-1
		return (j === 0 ? 0.5
			: j === M ? (N - 0.5)
			: rs[i] - 0.5 + (rs[j] - rs[i]) * (x - vs[i]) / (vs[j] - vs[i])
		) / N
	},
	/**
	 * @param {number} x - probability or array of probabilities
	 * @return {number} value or array of values
	 */
	f: function(x) {
		var vs = this.vs,
				rs = this.rs,
				M = Math.min(rs.length, rs[rs.length-1]),
				N = rs[M-1]
		if (x === vs[0] || x === vs[M-1]) return 0.5/N
		var j = topIndex(vs, x, M),
				i = j-1
		return j === 0 || j === M ? 0 : (rs[j] - rs[i]) / (vs[j] - vs[i]) / N
	},
	push: function(x) {
		var vs = this.vs,
				rs = this.rs,
				M = Math.min(rs.length, rs[rs.length-1]),
				j = topIndex(this.vs, x, M),
				ir = 0
		// lossless insert
		if (M < vs.length) {
			for (ir=M; ir>j; --ir) {
				rs[ir] = ir+1
				vs[ir] = vs[ir-1]
			}
			rs[j] = j ? rs[j-1] + 1 : 1
			vs[j] = x
			//console.log(j, rs)
			if (M!==rs.length-1) ++rs[rs.length-1]
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

function topIndex(arr, v, max) {
	var low = 0
	while (low < max) {
		var mid = (low + max) >>> 1
		if (arr[mid] < v) low = mid + 1
		else max = mid
	}
	return max
}
