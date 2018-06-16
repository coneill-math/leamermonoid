# leamermonoid
Provides a Python class `LeamerMonoid` for use with the computer algebra system [Sage](http://sagemath.org/) that adds functionality for working with Leamer monoids.  Internally, it uses the [GAP](http://www.gap-system.org/) package [numericalsgps](http://www.gap-system.org/Packages/numericalsgps.html) via the [NumericalSemigroup](https://github.com/coneill-math/numsgps-sage) class.  Several factorization invariants are supported, and more functionality will continue to be added.  

You can find action shots in the `images` folder.  

Please note that this is an *alpha version* and subject to change without notice.  

## License
leamermonoid is released under the terms of the [MIT license](https://tldrlegal.com/license/mit-license).  The MIT License is simple and easy to understand and it places almost no restrictions on what you can do with this software.

## Usage
To use this class, you must first install [numsgps-sage](https://github.com/coneill-math/numsgps-sage).  Next, simply place `LeamerMonoid.sage` in the same directory as `NumericalSemigroup.sage`.  

The following code fragment gives an overview of how to use the `LeamerMonoid` class from within Sage, and more complete documentation will be added in the near future.

	load('/PATH_TO_FILES/NumericalSemigroup.sage')
	load('/PATH_TO_FILES/LeamerMonoid.sage')
	S = LeamerMonoid([13,17,22,40],4)
	S.Plot().show()
	print S.frob
	print S.LengthSet(120,15)
	print S.DeltaSet(120,15)
	print S.CatenaryDegree(120,15)
