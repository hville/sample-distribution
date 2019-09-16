var testSample = require('./sample-tester')

module.exports = function(factories, randomFcns, N, M, S) {
	var keys = Object.keys(factories),
			fcns = keys.map(k => factories[k]),
			ress = keys.map(k => ({key:k, err:0, rms:0, ops: 0}))

	randomFcns.forEach(rfn => {
		for (var si=0; si<S; ++si) {
			for (var i=0, samples=Array(N); i<N; ++i) samples[i] = rfn(i)
			testSample(samples, fcns.map(f => f(M)), ress)
		}
	})
	ress.forEach(res => {
		res.ops = Math.round(1e2*N*S/res.ops*randomFcns.length)
		res.err = Math.round(1e6*res.err/N/S/randomFcns.length)
		res.rms = Math.round(1e4*Math.sqrt(res.rms/N/S/randomFcns.length))
	})
	return ress
}
