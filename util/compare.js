var comp = require('./comparator'),
		REC = require('..'),
		HDG = require('h-digest'),
		TDG = require('tdigest'),
		iZ = require('norm-dist/icdf')

var N = 2000,
		M = 20,
		S = 200

var randomFcns = [
	i => Math.exp(iZ(Math.random())+0.2), //logp
	i => -Math.exp(iZ(Math.random())+0.2),//logn
	i => iZ(Math.random()), // norm
	i => Math.floor(iZ(Math.random())*50), //disc
	i => iZ(Math.random()) + Math.random() < 0.3 ? -2 : +2, //bi-mode
	i => Math.random()*i, //min-max
	i => i, //sort
	i => -i //rvrs
]

var recs = {
	keepAve: () => new REC(M),
	hdigest: () => {
		var rec = HDG(M)
		rec.R = function(r) { return this.quantile( r/(N+1) ) }
		return rec
	},
	tdigest: () => {
		var rec = new TDG.TDigest(0.8, M, 1.1)
		rec.R = function(r) { return this.percentile( r/(N+1) ) }
		return rec
	}
}

function pad(ress, len) {
	return ress.reduce(
		(msg, res) => msg + `${res.key.padEnd(9)}  err: ${res.err.toFixed(0).padStart(len)}  rms: ${res.rms.toFixed(0).padStart(len)}  ops: ${res.ops.toFixed(0).padStart(len)}\n`,
		''
	)
}

var ress = comp(recs, randomFcns, N, M, S)
console.log(pad(ress, 8))
console.log('err: ', ress.sort( (a,b) => Math.abs(a.err) - Math.abs(b.err) ).map(v=>v.key).join(' < '))
console.log('rms: ', ress.sort( (a,b) => a.rms - b.rms ).map(v=>v.key).join(' < '))
console.log('ops: ', ress.sort( (a,b) => b.ops - a.ops ).map(v=>v.key).join(' < '))
