# polynomial.js
JavaScript file for working with polynomials over finite fields

License is omitted for the time being, until someone finds the files useful.

The use of the term "infinite fields" for {0, 1, 2, ...} is just a matter of programming convenience. It is mathematically *wrong* because there does not exist an additive inverse or a multiplicative inverse for some elements.

**polynomial.js** is the main file. **tests.html** checks for error using the console. No error should be returned in the console if any part of **polynomial.js** were to be changed.

Copyright 2016 @litena

# Documentation
Polynomials are created using function declaration, not via `new`.

* `polynomial(inputArray, inputIndeterminate, inputGalois)`
  * Create a polynomial object. `inputArray` hold the coefficients of the polynomial, with the value at index *n* being the coefficient of the term with exponent *n*. `inputIndeterminate` takes a string describing the indeterminate of the polynomial. `inputGalois` denotes the characteristic of the field.

# Available functions.

## `.array()`
### Returns the array of coefficients, with the index *n* being the coefficient of the term with exponent *n*. Only integers are allowed. Values are automatically changed to be between {0, 1, ..., *p*}, where *p* is the characteristic of the field.
