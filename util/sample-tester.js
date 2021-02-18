export default function(samples, recs, results) {
	var ress = results || recs.map(() => ({err:0, rms:0, ops: 0})),
			vals = normalize(samples)
	recs.forEach( (rec,i) => {
		var hr = process.hrtime()
		vals.forEach(rec.push, rec)
		hr = process.hrtime(hr)
		ress[i].ops += hr[0]/1e6 + hr[1]*1000
	})
	vals.sort((a, b) => a-b)
	recs.forEach((rec, ri) => {
		var res = ress[ri],
				err = 0
		for (var i=0; i<vals.length; ++i) {
			err = rec.R(i+1) - vals[i]
			res.err += err
			res.rms += err * err
		}
	})
	return ress
}

function normalize(samples) {
	var s1 = 0,
			s2 = 0,
			N = samples.length
	for (let v of samples) {
		s1 += v
		s2 += v*v
	}
	s2 = Math.sqrt( ( s2 - s1*s1/N ) / (N-1) ) //sample variance
	s1 /= N
	return samples.map(v => (v-s1) / s2)
}
