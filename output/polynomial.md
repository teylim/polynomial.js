# Global





* * *

### polynomial(inputArray, inputIndeterminate, inputGalois) 

**Parameters**

**inputArray**: `Array`, the integer coefficients of the polynomial. The value at index *n* is the coefficient of the term with exponent *n*. In infinite fields, the coefficients must be non-negative.

**inputIndeterminate**: `String`, the indeterminate of the polynomial. The default value is *x*.

**inputGalois**: `Number`, the prime (or *Infinity*) characteristic of the field. The default value is *Infinity*. When *Infinity* is specified, the coefficients take values in *{0, 1, 2, ...}*.

**Returns**: `polynomial`, the polynomial.


## Class: polynomial


### polynomial.array() 

**Returns**: `Array`, the coefficients of the polynomial, with the index *n* being the coefficient of the term with exponent *n*. Values are automatically changed to be between {0, 1, ..., *p*}, where *p* is the characteristic of the field.

### polynomial.indeterminate() 

**Returns**: `String`, the indeterminate of the polynomial.

### polynomial.galois() 

**Returns**: `Number`, the characteristic of the field.

### polynomial.clone() 

**Returns**: `polynomial`, a copy of the polynomial.

### polynomial.reciprocal() 

**Returns**: `polynomial`, the reciprocal polynomial.

### polynomial.degree() 

**Returns**: `Number`, the degree of the polynomial.

### polynomial.coefficient(n) 

**Parameters**

**n**: `Number`, the index.

**Returns**: `Number`, the coefficient of the term with exponent `n`.

### polynomial.leading() 

**Returns**: `Number`, the leading coefficient of the polynomial.

### polynomial.isMonic() 

**Returns**: `Boolean`, the polynomial is monic.

### polynomial.isZero() 

**Returns**: `Boolean`, the polynomial is *0*.

### polynomial.print() 

**Returns**: `String`, a formatted string of the polynomial, e.g., `GF(2): x ^ 2 + x + 1`.

### polynomial.differentiate() 

**Returns**: `polynomial`, the formal derivative of the polynomial.

### polynomial.isBinaryOperable(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `Object`, the value with key `error` shows if the polynomials have different indeterminates or different orders of fields.

### polynomial.add(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial

**Returns**: `polynomial`, the polynomial obtained when `inputPolynomial2` is added to the polynomial, if possible. Otherwise, returns an object with key `error`.

### polynomial.subtract(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `polynomial`, the polynomial obtained when `inputPolynomial2` is subtracted from the polynomial, if possible. Otherwise, returns an object with key `error`.

### polynomial.multiply(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `polynomial`, the polynomial obtained when `inputPolynomial2` is multiplied to the polynomial, if possible. Otherwise, returns an object with key `error`.

### polynomial.power(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `polynomial`, the product of `n` repetitions of the polynomial, if possible. Otherwise, returns an object with key `error`.

### polynomial.givesQuotient(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `polynomial`, the *largest* polynomial which product with `inputPolynomial2` is *smaller than* or *equal to* the polynomial, if possible. Otherwise, returns an object with key `error`.

### polynomial.givesRemainder(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `polynomial`, the *smallest* non-negative polynomial obtained when `inputPolynomial2` and its products with powers of *x* are repeatedly subtracted from the polynomial, if possible. Otherwise, returns an object with key `error`.

### polynomial.of(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `polynomial`, the polynomial obtained when `inputPolynomial2` is substituted at every instance of the indeterminate. The substituted indeterminate is evaluated in the original field.

### polynomial.additiveOrderOf(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `Number`, the smallest positive integer *k* where the `inputPolynomial2` multiplied by *k* gives *0*, if possible.

### polynomial.multiplicativeOrderOf(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `Number`, the smallest positive integer *k* where `inputPolynomial2` raised to the power of *k* gives *1*, if possible.

### polynomial.hasPrimitiveElement(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `Boolean`, `inputPolynomial2` is a primitive element of the polynomial

### polynomial.hasFactor(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `Boolean`, `inputPolynomial2` is a factor of the polynomial

### polynomial.hasRoot(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `Boolean`, the polynomial is *0* when evaluated at `inputPolynomial2`.

### polynomial.successor() 

**Returns**: `polynomial`, the smallest polynomial which is *greater than* the polynomial.

### polynomial.isLessThan(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `Boolean`, `inputPolynomial2` can be obtained by finding the recursive successor of the polynomial.

### polynomial.isGreaterThan(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `Boolean`, the polynomial can be obtained by finding the recursive successor of `inputPolynomial2`.

### polynomial.isEqualTo(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `Boolean`, the polynomial is neither "less than" nor "greater than" `inputPolynomial2`.

### polynomial.additiveGroupOf(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `Array`, the additive group of `inputPolynomial2` modulo the polynomial.

### polynomial.multiplicativeGroupOf(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `Array`, the multiplicative group of `inputPolynomial2` modulo the polynomial.

### polynomial.primitiveElements() 

**Returns**: `Array`, the primitive elements of the polynomial.

### polynomial.factors() 

**Returns**: `Array`, the factors of the polynomial. The polynoial 0 has infinitely many factors.* For infinite fields, returns  * "undefined" for non-constant polynomials,  * prime factorization for constant polynomials larger than *1*,  * "undefined" for *1* and  * "undefined" for *0*.* For finite fields, returns  * irreducible polynomials factorization for constant polynomials with positive degree,  * "undefined" for constant polynomials larger than *1*,  * empty array for *1* and   * "undefined" for *0*.

### polynomial.isIrreducible() 

**Returns**: `Boolean`, the polynomial is irreducible (but possibly non-monic).

### polynomial.isMinimalPolynomialOf(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `Boolean`, the polynomial in a field is the minimal polynomial of `inputPolynomial2` in the extension field.

### polynomial.primitivePolynomialsInExtensionOf(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `Array`, the primitive polynomials corresponding to the primitive elements returned by `.primitiveElements()`.

### polynomial.standard(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `Array`, the serialisation counting the occurences of a factor starting with the *0* followed recursively by successors.

### polynomial.fromStandard(standardArray) 

**Parameters**

**standardArray**: `Array`, a standard array.

**Returns**: `polynomial`, the deserialisation of `standardArray`, taking `indeterminate()` and `galois()` from the polynomial.

### polynomial.gcd(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `polynomial`, the greatest common divisor of the polynomial and `inputPolynomial2`.

### polynomial.lcm(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `polynomial`, the lowest common multiple of the polynomial and `inputPolynomial2`.

### polynomial.cycles(inputPolynomial2) 

**Parameters**

**inputPolynomial2**: `polynomial`, another polynomial.

**Returns**: `Array`, the polynomials *p, xp, xxp, ...*, where *p* is `inputPolynomial2` and *x* is once the indeterminate of the polynomial. The degree of the last element is *1* less than the degree of the polynomial.

### polynomial.cyclicCodeBasis(codeLength, codeDimension) 

**Parameters**

**codeLength**: `Number`, the length of the desired code.

**codeDimension**: `Number`, the dimension of the desired code.

**Returns**: `Array`, the arrays corresponding to *a* basis of a cyclic code with code length `codeLength` and code dimension `codeDimension`, taking galois() from the polynomial.



* * *










