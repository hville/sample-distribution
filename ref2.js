export default class D {

	constructor(size=32) {
		const vs = size.buffer ? size : new Float64Array(size.byteLength ? size : 2*size)
		// make properties !writeable !configurable and !enumerable
		Object.defineProperties(this, {
			vs: {value: vs},
			rs: {value: new Float64Array(vs.buffer, vs.byteOffset + vs.byteLength/2, vs.length/2)}
		})
	}

	// Number of samples
	get N() { return this.rs[this.rs.length-1] }

	// Expected Value
	get E() { return this.Σ(1) / this.N }

	// Sample Variance
	get V() {
		const N = this.N
		return ( this.Σ(2) - this.Σ(1)**2 / N ) / (N-1)
	}

	// Sample Standard Deviation
	get S() {
		const v = this.V
		return v < 0 ? 0 : Math.sqrt(v)
	}

	/**
	 * Σ(X**p)
	 * exact when there is no compression
	 * with compression, range between values treated as a uniform distribution
	 *
	 * @param {number} order
	 * @return {number} Σ( X^pow )
	 */
	Σ(pow) { //values as-is with internal uniform interval
		const vs = this.vs,
					rs = this.rs,
					N = Math.min(rs.length, rs[rs.length-1]), // in case the buffer is not full
					Mm = N-1,
					Op = pow + 1
		if (pow === 0) return rs[Mm]
		if (pow === 1) { //same as below but simplified
			let sum = vs[0] + vs[Mm] // correction at edge to match actual discrete result (PAD cancels out)
			for (let i=0; i<Mm; ++i) sum += (vs[i+1] + vs[i]) * (rs[i+1] - rs[i])
			return sum / Op
		}
		let sum = vs[0]**pow
		for (let i=1; i<N; ++i) {
			// https://en.wikipedia.org/wiki/Continuous_uniform_distribution#Moments
			sum += vs[i]**pow
					+ (rs[i] - rs[i-1] - 1) * (vs[i]**Op - vs[i-1]**Op) / (vs[i] - vs[i-1]) / Op
		}
		return sum
	}

	/**
	 * Origin Moments
	 * https://en.wikipedia.org/wiki/Continuous_uniform_distribution#Moments
	 *
	 * @param {number} order
	 * @return {number} E( X^order )
	 */
	M(order) { return this.Σ(order) / this.N }

	/**
	 * Quantile function, provide the value for a given probability
	 * @param {number} prob - probability or array of probabilities
	 * @return {number} value or array of values
	 */
	Q(prob) {
		const vs = this.vs,
					rs = this.rs,
					M = Math.min(rs.length, rs[rs.length-1]), // in case the buffer is not full
					h = rs[M-1] * prob + 0.5, // 0.5 <= h <= N+0.5
					j = topIndex(rs, h, M), //      0 <= j <= M
					i = j-1
		return j === 0 ? vs[0]
			: j === M ? vs[M-1]
			: vs[i] + (vs[j] - vs[i]) * (h-rs[i]) / (rs[j]-rs[i])
	}

	/**
	 * @param {number} x - probability or array of probabilities
	 * @return {number} value or array of values
	 */
	F(x) {
		const vs = this.vs,
					rs = this.rs,
					M = Math.min(rs.length, rs[rs.length-1]), // in case the buffer is not full
					N = rs[M-1],
					j = topIndex(vs, x, M),
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
	f(x) {
		const vs = this.vs,
					rs = this.rs,
					M = Math.min(rs.length, rs[rs.length-1]),
					N = rs[M-1]
		if (x === vs[0] || x === vs[M-1]) return 0.5/N
		const j = topIndex(vs, x, M),
					i = j-1
		return j === 0 || j === M ? 0 : (rs[j] - rs[i]) / (vs[j] - vs[i]) / N
	}

	/**
	 * @param {Object} ctx - canvas 2D context
	 * @param {number} xMin
	 * @param {number} xMax
	 * @return {void}
	 */
	plotF(ctx, xMin = this.vs[0], xMax = this.vs[this.rs.length-1])	{
		const rs = this.rs,
					vs = this.vs,
					xScale = ctx.canvas.width / (xMax-xMin),
					yScale = ctx.canvas.height / (rs[rs.length-1]),
					H = ctx.canvas.height

		ctx.beginPath()
		ctx.moveTo( (vs[0]-xMin)*xScale, H )
		for (let i=0; i<rs.length; ++i) ctx.lineTo( (vs[i]-xMin)*xScale, H-(rs[i]-0.5)*yScale )
		ctx.lineTo( (vs[rs.length-1]-xMin)*xScale, 0 )
	}
	/**
	 * @param {Object} ctx - canvas 2D context
	 * @param {number} xMin
	 * @param {number} xMax
	 * @param {number} yMax
	 * @return {void}
	 */
	plotf(ctx, xMin = this.vs[0], xMax = this.vs[this.rs.length-1], yMax = 5/(this.vs[this.rs.length-2]-this.vs[1])) {
		const rs = this.rs,
					vs = this.vs,
					xScale = ctx.canvas.width / (xMax-xMin),
					yScale = ctx.canvas.height * (this.vs[this.rs.length-1]-this.vs[0]) / (6*rs[rs.length-1]),
					H = ctx.canvas.height

		let x = (vs[0]-xMin) * xScale
		ctx.beginPath()
		ctx.moveTo( x, H )
		for (let i=0, j=1, y=0; j<rs.length; i=j++) {
			ctx.lineTo( x, y = H - (rs[j]-rs[i])/(vs[j]-vs[i]) * yScale ) //up
			ctx.lineTo( x = (vs[j]-xMin) * xScale, y ) //right
		}
		ctx.lineTo( x, H )
	}

	/**
	 * Adds a value, compressed only if buffer full
	 * @param {number} x
	 */
	push(x) {
		const vs = this.vs,
					rs = this.rs,
					M = Math.min(rs.length, rs[rs.length-1])
		let j = topIndex(this.vs, x, M)
		// lossless insert
		if (M < rs.length) {
			for (let ir=M; ir>j; --ir) {
				rs[ir] = ir+1
				vs[ir] = vs[ir-1]
			}
			rs[j] = j ? rs[j-1] + 1 : 1
			vs[j] = x
			if (M!==rs.length-1) ++rs[rs.length-1]
		}
		// compression, droping a value while maintaining the interval average
		else if (j === M) {
			--j
			const i = j-1,
						h = i-1,
						Δwv = vs[j]-vs[i],
						Δxu = x - vs[h],
						rjh = rs[i]*(vs[j]-vs[h])
			if (Δxu !== 0) {
				const r_w = ( rs[j]*(x-vs[i]) + rjh - rs[h]*Δwv ) / Δxu,
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
		else if (j === 0) {
			const u = vs[0],
						Δwx = vs[2] - x
			if (Δwx !== 0) {
				//DROP u = vs[0] : r'(w-x) = (w+v-2x) + r(w-u)          //KEEP v
				//KEEP u = vs[0] : r'(w-x) = (w+v-2x) + r(w-u) - s(v-u) //DROP v
				const r_v = (vs[2] + vs[1] - 2*x + rs[1]*(vs[2] - u)) / Δwx, // DROP u KEEP v
							r_u = r_v - rs[2] * (vs[1]-u) / Δwx                    // KEEP u DROP v
				if ( (u + vs[1] > x + vs[2] && r_u > rs[0]) || r_v > rs[2]+1) {
					vs[1] = u
					rs[1] = r_u
				}
				else rs[1] = r_v
				vs[0] = x
			}
			for (let ir=2; ir<rs.length; ++ir) ++rs[ir]
		}
		else if ( j !== 1 && (j === M-1 || 2*x < vs[j]+vs[j-1] ) ) {
			// u < v < x < w
			--j
			let k = j+1,
					i = j-1
			const w = vs[k],
						v = vs[j],
						Δwu = w-vs[i]
			if (Δwu !== 0) {
				const r_x = rs[j] + ( w-x + (x-v)*(rs[k] - rs[i]) ) / Δwu
				if ( vs[i]+w < v+x || r_x > rs[k]+1) rs[j] += ( w+v-2*x ) / Δwu
				else {
					rs[j] = r_x
					vs[j] = x
				}
			}
			while(++j<rs.length) ++rs[j]
		}
		else {
			// u < x < v < w
			let k = j+1,
					i = j-1
			const w = vs[k],
						v = vs[j],
						Δwu = w-vs[i]
			if (Δwu !== 0) {
				const r_x = rs[j] + ( w-x + (x-v)*(rs[k] - rs[i]) ) / Δwu
				if ( v+x < vs[i] + w || r_x < rs[i]) rs[j] += ( w+v-2*x ) / Δwu
				else {
					rs[j] = r_x
					vs[j] = x
				}
			}
			while(++j<rs.length) ++rs[j]
		}
	}
}

function topIndex(arr, v, max) {
	let low = 0
	while (low < max) {
		const mid = (low + max) >>> 1
		if (arr[mid] < v) low = mid + 1
		else max = mid
	}
	return max
}

function getSum2(vs,rs,i,j,k) {
	return ( vs[k]+vs[j] )*( rs[k]-rs[j] ) + ( vs[j]+vs[i] )*( rs[j]-rs[i] )
}
function rjFix(vs,rs,i,j,k,sum2) {
	rs[j] = ( (vs[k]+vs[j])*rs[k] - (vs[j]+vs[i])*rs[i] - sum2 ) / (vs[k]-vs[i])
	//if (Math.abs(sum2-getSum2(vs,rs,i,j,k)) > 0.001) throw Error
}
