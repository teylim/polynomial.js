# polynomial.js
JavaScript file for working with polynomials over finite fields

License is omitted for the time being, until someone finds the files useful.

The use of the term "infinite fields" for {0, 1, 2, ...} is just a matter of programming convenience. It is mathematically *wrong* because there does not exist an additive inverse or a multiplicative inverse for some elements.

**polynomial.js** is the main file. **tests.html** checks for error using the console. No error should be returned in the console if any part of **polynomial.js** were to be changed.

Copyright 2016 @litena

# Documentation
Polynomials are created using function declaration, not via `new`.

* `polynomial(inputArray, inputIndeterminate, inputGalois)`
  * Create a polynomial object. `inputArray` hold the integer coefficients of the polynomial, with the value at index *n* being the coefficient of the term with exponent *n*. `inputIndeterminate` takes a string describing the indeterminate of the polynomial. `inputGalois` denotes the prime (or Infinity\*) characteristic of the field.
\* When *Infinity* is specified, the coefficients take values in {0, 1, 2, ...}
# Available functions.

* `.array()`
  * Returns the array of coefficients, with the index *n* being the coefficient of the term with exponent *n*. Values are automatically changed to be between {0, 1, ..., *p*}, where *p* is the characteristic of the field.
* `.indeterminate()`
  * Returns the indeterminate of the polynomial. The default value is "x".
* `.galois()`
  * Returns the characteristic of the field.
* `.clone()`
  * Returns a duplicate of the polynomial.
* `.reciprocal()`
  * Returns the reciprocal polynomial.
* `.degree()`
  * Returns the degree of the polynomial.
* `.coefficient(n)`
  * Returns the coefficient of the term with exponent `n`.
* `.leading()`
  * Returns the leading coefficient the polynomial.
* `.isMonic()`
  * Returns `true` if the polynomial is monic and `false` otherwise.
* `.isZero()`
  * Returns `true` if the polynomial is the zero polynomial and `false` otherwise.
* `.print()`
  * Obtained a formatted string of the polynomial, e.g., `GF(2): x ^ 2 + x + 1`.
* `.differentiate()`
  * Returns the formal derivative of the polynomial.
* `.isBinaryOperable(inputPolynomial2)`
  * Returns an object with key `error` if the polynomials have the different indeterminates or different order of field.
* `.add(inputPolynomial2)`
  * Returns the polynomial obtained when `inputPolynomial2` is added to the polynomial. Otherwise, returns an object with key `error`.
* `.subtract(inputPolynomial2)`
  * Returns the polynomial obtained when `inputPolynomial2` is subtracted from the polynomial. Otherwise, returns an object with key `error`.
* `.multiply(inputPolynomial2)`
  * Returns the polynomial obtained when `inputPolynomial2` is multiplied to the polynomial. Otherwise, returns an object with key `error`.
* `.power(n)`
  * Returns the polynomial obtained when polynomial is product of `n` repetitions of the polynomial. Otherwise, returns an object with key `error`.
* `.givesQuotient(inputPolynomial2)`
  * Returns the largest polynomial which product with `inputPolynomial2` is smaller than the polynomial. Otherwise, returns an object with key `error`.
* `.givesRemainder(inputPolynomial2)`
  * Returns the smallest non-negative polynomial obtained when `inputPolynomial2` is repeatedly subtracted from the polynomial. Otherwise, returns an object with key `error`.
* `.of(inputPolynomial2)`
  * Returns the polynomial obtained when `inputPolynomial2` is substituted at every instance of the indeterminate. The substituted indeterminate is evaluated in the original field.
* `.orderOf(inputPolynomial2)`
  * Returns the smallest *k* where `inputPolynomial2` raised to the power of *k* gives 1, if possible.
* `.hasPrimitiveElement(inputPolynomial2)`
  * Returns `true` if `inputPolynomial2` is a primitive element of the polynomial and `false` otherwise.
* `.hasFactor(inputPolynomial2)`
  * Returns `true` if `inputPolynomial2` is a factor of the polynomial and `false` otherwise.
* `.hasRoot(inputPolynomial2)`
  * Returns `true` if the polynomial is zero when evaluated at `inputPolynomial2` and `false` otherwise.
 
