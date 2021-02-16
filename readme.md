<!-- markdownlint-disable MD004 MD007 MD010 MD041 MD022 MD024 MD032 -->
# sample distribution

*takes a large and continuous data stream and continuously only retain a reduced [empirical CDF](https://en.wikipedia.org/wiki/Empirical_distribution_function) approximation*

***small, simple, no dependencies***

• [Example](#example) • [Features](#features) • [Limitations](#limitations) • [Why](#why) • [API](#api) • [License](#license)

# Example

```javascript
import Recorder from 'sample-distribution'

// recorder with 5 retained samples
var recorder = new Recorder(5)
[5,3,6,2,7,1,8].forEach(recorder.push, recorder)

console.log('minimum:', recorder.Q(0)) // minimum:0
console.log('median:', recorder.Q(0.5)) // median:4
console.log('maximum:', recorder.Q(1)) // maximum:8
```

# Features

* constant memory use, no compression steps and/or triggered garbage collection
* significantly faster than other implementations (about 3-5x faster)
* mean-preserving compression (empirical cdf with same average as samples)
* No value interpolation. Values are kept as-is, but ranks are interpolated.

# Limitations

* other than the mean, the other moments (variance, skew, kurtosis) are not preserved
* The moments are those of the generated curve that approximate those of the source samples

# Background and related projects

* [tdigest](https://www.npmjs.com/package/tdigest) based on the [work of Dunning](https://github.com/tdunning/t-digest)
* [h-digest](https://www.npmjs.com/package/h-digest) based on the above but keeping values and adjusting ranks only

This module attempts to match the underlying CDF as closely as possible by adjusting local ranks to keep the local average when discarding values.
By some measure, the root mean square error for each original sample value is 10 times better than the 2 other implementations above.
Type `npm run compare` or see `./util/compare.js` for a benchmark and error comparison for multiple distribution types.

# API

## Creation

* `var recorder = new Recorder(L)`: creates a recorder that will keep `L` values and associated `L` ranks

## Properties - read-only getters
* `.N` number: total samples received
* `.E` number: average of samples received and of the resulting approximated cdf
* `.S` number: standard variation of the resulting aproximated cdf
* `.V` number: variance of the resulting aproximated cdf

## Properties - internal
* `.vs` array: internal store of retained sample values
* `.rs` array: internal store of retained value approximated ranks

## Methods
* `.push(number)` void: sample value(s) to be added
* `.F(value:number)` number: cummulative probability (cdf) of the specified value
* `.f(value:number)` number:  probability density (pdf) of the specified value
* `.Q(probability:number)` number: estimated value for specified probability
* `.M(order:number)` number: estimated origin moment (ie. E = M(1), V=M(2)-E^2)

# License

[MIT](http://www.opensource.org/licenses/MIT) © [Hugo Villeneuve](https://github.com/hville)
