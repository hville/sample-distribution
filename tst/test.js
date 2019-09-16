var c = require('cotest'),
		tester = require('../util/sample-tester'),
		Rec = require('../')

c('len = 5, some identical values, non compressed', t => {
	var hd = new Rec(5)
	hd.push(1)
	hd.push(2)
	hd.push(3)
	t('==', hd.Q(0.5), 2)

	hd.push(2)
	t('==', hd.Q(0), 1)
	t('==', hd.Q(0.5), 2)
	t('==', hd.Q(1), 3)
	hd.push(2)
	t('==', hd.Q(0), 1)
	t('==', hd.Q(0.5), 2)
	t('==', hd.Q(1), 3)

})
c('len = 5, compressed, sorted (all new max)', t => {
	var hd = new Rec(5)
	testSet(t, hd, [0,1,2,3,4,5,6,7,8,9])
})
c('len = 5, reverse sorted (all new min)', t => {
	var hd = new Rec(5)
	testSet(t, hd, [8,7,6,5,4,3,2,1,0])
})
c('len = 5, mix-min-max', t => {
	testSet(t, new Rec(5), [4,5,3,6,2,7,1,8,0])
})
c('len = 9, 101 sorted', t => {
	var n = 99,
			rnd = []
	for (var i=0; i<n; ++i) rnd.push(i)
	testSet(t, new Rec(9), rnd)
})
c('len = 9, 101 reversed', t => {
	var rnd = []
	for (var i=0; i<101; ++i) rnd.push(101-i-1)
	testSet(t, new Rec(9), rnd)
})
c('uniform random, m=9, n=1000', t => {
	var	n = 1001,
			rnd = [],
			hd = new Rec(9)
	for (var i=0; i<n; ++i) rnd.push(Math.random())
	testSet(t, hd, rnd)
})
c('skewed random, m=9, n=1000', t => {
	var rnd = []
	for (var i=0; i<1001; ++i) rnd.push(Math.random()+Math.random()*Math.random())
	testSet(t, new Rec(9), rnd)
})
c('lots of repeated values n=9, n=1000', t => {
	var rnd = []
	for (var i=0; i<1000; ++i) rnd.push(Math.random() < 0.2 ? 0 : Math.random()*Math.random())
	testSet(t, new Rec(9), rnd)
})
c('RMS', t => {
	var	trd = 333,
			n = 3*trd,
			rnd = Array(n),
			rms = [0,0,0,0],
			sum = (r,v) => r+v
	function inc(a, b) {return a-b}
	for (var i=0; i<n; ++i) {
		var hd = new Rec(21)
		for (var j=0; j<n; ++j) {
			rnd[j] = Math.random()-0.5
			hd.push(rnd[j])
		}
		rnd.sort(inc)
		var k = 0
		while (k < trd) rms[0] += Math.pow(hd.Q((k-0.5)/n) - rnd[k++], 2)
		while (k < 2*trd) rms[1] += Math.pow(hd.Q((k-0.5)/n) - rnd[k++], 2)
		while (k < n) rms[2] += Math.pow(hd.Q((k-0.5)/n) - rnd[k++], 2)
		rms[3] += Math.pow(hd.E * n - rnd.reduce(sum), 2)
	}
	t('<', rms[0]/n/trd, 0.001, 'low rms')
	t('<', rms[1]/n/trd, 0.001, 'mid rms')
	t('<', rms[2]/n/trd, 0.001, 'top rms')
	t('<', rms[3]/n, 1e-12, 'ave'+rnd.reduce(sum)/n+' '+hd.ave)
})
function testSet(t, cdf, set) {
	var results = tester(set,[cdf])[0],
			N = set.length
	results.err /= N
	results.rms = Math.sqrt(results.rms/N)
	t('<', Math.abs(cdf.E), 1e-12, 'mean')
	t('<', Math.abs(results.err), 2e-5, 'bias')
	t('<', Math.abs(results.rms), 2e-1, 'rms')
	t('<', Math.abs(cdf.F(cdf.Q(0.5))-0.5), 1e-6, 'F(Q(0.5)) ~= 0.5')
	t('<', Math.abs(cdf.Q(cdf.F(0))), 1e-6, 'Q(F(0)) ~= 0')
	t('===', cdf.R(0), cdf.Q(0), 'R(0) === Q(0)')
	t('===', cdf.R(N), cdf.Q(1), 'R(N) === Q(1)')
}
