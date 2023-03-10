export default class D {

	constructor(size=32) {
		const buffer = size.buffer || (size.byteLength ? size : new ArrayBuffer(size << 4)),
					offset = size.byteOffset || 0,
					byteLn = (size.byteLength || buffer.byteLength) >> 1,
					length = byteLn >> 3
		Object.defineProperties(this, {
			vs: {value: new Float64Array(buffer, offset, length)},
			rs: {value: new Float64Array(buffer, offset + byteLn, length)}
		})
	}
	// for transfers and copies
	get data() { return new DataView(this.vs.buffer, this.vs.byteOffset, this.vs.byteLength << 1) }

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
					M = Math.min(rs.length, rs[rs.length-1]), // in case the buffer is not full
					Mm = M-1,
					Op = pow + 1
		if (pow === 0) return rs[Mm]
		if (pow === 1) { //same as below but simplified
			let sum = vs[0] + vs[Mm] // correction at edge to match actual discrete result (PAD cancels out)
			for (let i=0; i<Mm; ++i) sum += (vs[i+1] + vs[i]) * (rs[i+1] - rs[i])
			return sum / Op
		}
		let sum = vs[0]**pow
		for (let i=1; i<M; ++i) {
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
		const j = topIndex(vs, x, M)
		return j === 0 || j === M ? 0 : (rs[j] - rs[j-1]) / (vs[j] - vs[j-1]) / N
	}

	/**
	 * @param {Object} ctx - canvas 2D context
	 * @param {number} vMin
	 * @param {number} vMax
	 * @return {void}
	 */
	plotF(ctx, vMin=this.vs[0], vMax=this.vs[this.rs.length-1])	{
		const rs = this.rs,
					vs = this.vs,
					xScale = (ctx.canvas.width-1) / (vMax-vMin),
					yScale = (ctx.canvas.height-1) / (rs[rs.length-1]),
					H = ctx.canvas.height,
					getX = v => 0.5 + Math.round((v-vMin) * xScale),
					getY = r => H-0.5 - Math.round(r * yScale)
		ctx.beginPath()
		ctx.moveTo( getX( Math.min(vs[0], vMin) ), H-0.5 )
		ctx.lineTo( getX( vs[0] ), H-0.5 )
		for (let i=0; i<rs.length; ++i) ctx.lineTo( getX( vs[i] ), getY( rs[i] ) )
		ctx.lineTo( getX( vs[rs.length-1] ), 0.5 )
		ctx.lineTo( getX( Math.max( vs[rs.length-1], vMax ) ), 0.5 )
	}
	/**
	 * @param {Object} ctx - canvas 2D context
	 * @param {number} vMin
	 * @param {number} vMax
	 * @param {number} yMax
	 * @return {void}
	 */
	plotf(ctx, vMin=this.vs[0], vMax=this.vs[this.rs.length-1], yMax = 1/(this.Q(0.75)-this.Q(0.25)) ) {
		// f(mode)*IQR ~= 0.5(uniform) 0.54(normal) 0.55(logistic) 0.59(triangular) 0.64(cauchy) 0.69(laplace)
		// aim for yMax ~= 3*E(f(mode)) ~= 3*0.5/IQR
		const rs = this.rs,
					vs = this.vs,
					xScale = (ctx.canvas.width-1) / (vMax-vMin),
					yScale = (ctx.canvas.height-1) / yMax / rs[rs.length-1],
					H = ctx.canvas.height,
					getX = v => 0.5 + Math.round((v-vMin) * xScale),
					getY = drdv => H-0.5 - Math.round(drdv * yScale)
		let x = getX(Math.min(vs[0], vMin)),
				y = H
		ctx.beginPath()
		ctx.moveTo( getX( Math.min(vs[0], vMin) ), H-0.5 )
		ctx.lineTo( x = getX(vs[0]), H-0.5 ) //right
		for (let i=0, j=1; j<rs.length; i=j++) {
			ctx.lineTo( x, y = getY( (rs[j]-rs[i])/(vs[j]-vs[i]) ) ) //up||down
			ctx.lineTo( x = getX( vs[j] ), y ) //right
		}
		ctx.lineTo( x, H-0.5 ) //down
		ctx.lineTo( getX( Math.max(vs[rs.length-1], vMax) ), H-0.5 ) //right
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
		else if (j === M) { // new maximum
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
		else if (j === 0) { // new minimum
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
		else if ( j !== 1 && (j === M-1 || 2*x < vs[j+1]+vs[j-2] ) ) { //vs[j]+vs[j-1] ) ) {
			// not max, not min : u < v < x < w
			--j
			let k = j+1,
					i = j-1
			const w = vs[k],
						v = vs[j],
						Δwu = w-vs[i]
			if (Δwu !== 0) {
				const r_xΔwu = rs[j]*Δwu + ( w-x + (x-v)*(rs[k] - rs[i]) )
				if ( vs[i]+w < v+x || r_xΔwu >= (rs[k]+1)*Δwu) rs[j] += ( w+v-2*x ) / Δwu
				else {
					rs[j] = r_xΔwu/Δwu
					vs[j] = x
				}
			}
			while(++j<rs.length) ++rs[j]
		}
		else {
			// not max, not min : u < x < v < w
			let k = j+1,
					i = j-1
			const w = vs[k],
						v = vs[j],
						Δwu = w-vs[i]
			if (Δwu !== 0) {
				const r_xΔwu = rs[j]*Δwu + ( w-x + (x-v)*(rs[k] - rs[i]) )
				if ( x+v < vs[i]+w || r_xΔwu <= rs[i]*Δwu) rs[j] += ( w+v-2*x ) / Δwu
				else {
					rs[j] = r_xΔwu/Δwu
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
