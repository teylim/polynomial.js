var polynomial;

/**
 * @param {Array} inputArray the integer coefficients of the polynomial. The value at index *n* is the coefficient of the term with exponent *n*. In infinite fields, the coefficients must be non-negative.
 * @param {String} inputIndeterminate the indeterminate of the polynomial. The default value is *x*.
 * @param {Number} inputGalois the prime (or *Infinity*) characteristic of the field. The default value is *Infinity*. When *Infinity* is specified, the coefficients take values in *{0, 1, 2, ...}*.
 * @returns {polynomial} the polynomial.
 */
// IMPORTANT NOTE
// * The use of the term "infinite fields" for *{0, 1, 2, ...}* is just a matter of programming convenience. It is mathematically *wrong* because there does not exist an additive inverse or a multiplicative inverse for some elements.
// * Define the increase in degree across polynomials as increasing and between polynomials having the same degree, the ordering *0, 1, ...* for each coefficient as increasing.

polynomial = function (inputArray, inputIndeterminate, inputGalois) {

// Return object with key "error" if errors thrown
try{

  var modulo, zeros, dividedBy, inputArray, inputIndeterminate, inputGalois, i, array, indeterminate, galois, clone, reciprocal, degree, coefficient, leading, isMonic, isZero, print, differentiate, isBinaryOperable, add, subtract, multiply, power, givesQuotient, givesRemainder, of, additiveOrderOf, multiplicativeOrderOf, hasPrimitiveElement, hasFactor, hasRoot, successor, isLessThan, isGreaterThan, isEqualTo, additiveGroupOf, multiplicativeGroupOf, primitiveElements, factors, isIrreducible, isMinimalPolynomialOf, primitivePolynomialsInExtensionOf, standard, fromStandard, gcd, lcm, cycles, cyclicCodeBasis;

  //// PRIVATE FUNCTIONS

  modulo = function (m, n) {
    if (!isFinite(n)) {
      return m;
    }
    return ((m % n) + n) % n;
  };

  zeros = function (d) {
    var outputArray, i;
    outputArray = [];
    for (i = 0; i <= d; i += 1) {
      outputArray[i] = 0;
    }
    return outputArray;
  };

  dividedBy = function (inputPolynomial2) {
    var dividend, outputArray, d, i, c, temporaryArray;
    if (typeof isBinaryOperable(inputPolynomial2).error !== "undefined") {
      return isBinaryOperable(inputPolynomial2);
    }
    // inputPolynomial2 cannot be zero
    if (inputPolynomial2.isZero()) {
      throw "Division by zero";
    }
    // Perform division as if by hand
    dividend = clone();
    outputArray = degree() - inputPolynomial2.degree() >= 0 ? zeros(degree() - inputPolynomial2.degree()) : [];
    notFactorizable: for (d = degree(); d >= inputPolynomial2.degree(); d -= 1) {
      c = dividend.coefficient(d) || 0;
      while (c === 0) {
        d -= 1;
        if (d < inputPolynomial2.degree()) {
          break notFactorizable;
        }
        c = dividend.coefficient(d) || 0;
      }
      if (isFinite(galois())) {
        i = 1;
        while (modulo(i * inputPolynomial2.leading(), galois()) !== c) {
          i += 1;
        }
      } else {
        i = inputPolynomial2.leading() > 0 ? Math.floor(c / inputPolynomial2.leading()) : Math.ceil(c / inputPolynomial2.leading());
      }
      outputArray[d - inputPolynomial2.degree()] = i;
      // temporaryArray is an array of the same length as outputArray with only the leading entry not replaced by 0
      temporaryArray = zeros(d - inputPolynomial2.degree());
      temporaryArray[d - inputPolynomial2.degree()] = i;
      dividend = dividend.subtract(inputPolynomial2.multiply(polynomial(temporaryArray, indeterminate(), galois())));
    }
    return [polynomial(outputArray, indeterminate(), galois()), dividend];
  };

  //// MODIFY THE INPUT ARGUMENTS ////

  // Index 0 has exponent 0, index i has exponent 1, etc.
  // Scalar input is treated as constant polynomial
  if (typeof inputArray === "number") {
    inputArray = [inputArray];
  }
  // Default polynomial is 0
  if (!Array.isArray(inputArray)) {
    inputArray = [0];
  }
  // Coefficients must be integers
  for (i = 0; i < inputArray.length; i += 1) {
    if (inputArray[i] !== parseInt(inputArray[i], 10)){
      throw "Some coefficients are not integers";
    }
  }

  // Default indeterminate is "x"
  if (inputIndeterminate === "" || typeof inputIndeterminate !== "string") {
    inputIndeterminate = "x";
  }

  // Default order of field is infinite
  if (inputGalois <= 1 || typeof inputGalois !== "number") {
    inputGalois = Infinity;
  }
  // Order of finite field must be prime
  if (isFinite(inputGalois)) {
    if (inputGalois !== parseInt(inputGalois, 10)) {
      throw "Order of field is not an integer";
    }
    for (i = 2; i <= Math.sqrt(inputGalois); i += 1) {
      if (inputGalois % i === 0 && i !== inputGalois) {
        throw "Order of finite field is not prime";
      }
    }
  }
  // All integer coefficients of polynomials over infinite field must be nonegative
  if (!isFinite(inputGalois)) {
    for (i = 0; i < inputArray.length; i += 1) {
      if (inputArray[i] < 0){
        throw "Some coefficients are negative for finite field";
      }
    }
  }
  // Coefficients of polynomials over finite field must be integers within the finite field
  for (i = 0; i < inputArray.length; i += 1) {
    inputArray[i] = modulo(inputArray[i], inputGalois);
  }
  // Leading coefficient cannot be zero unless polynomial is 0
  while (inputArray[inputArray.length - 1] === 0) {
    inputArray.pop();
  }
  if (inputArray.length === 0) {
    inputArray = [0];
  }

  //// NO MORE CHANGES ARE TO BE MADE TO inputArray, inputIndeterminate AND inputGalois FROM HERE ON ////

  /**
   * @class polynomial
   */

  /**
   * @returns {Array} the coefficients of the polynomial, with the index *n* being the coefficient of the term with exponent *n*. Values are automatically changed to be between {0, 1, ..., *p*}, where *p* is the characteristic of the field.
   */
  array = function () {
    return inputArray.slice(0);
  };

  /**
   * @returns {String} the indeterminate of the polynomial.
   */
  indeterminate = function () {
    return inputIndeterminate;
  };

  /**
   * @returns {Number} the characteristic of the field.
   */
  galois = function () {
    return inputGalois;
  };

  /**
   * @returns {polynomial} a copy of the polynomial.
   */
  clone = function () {
    return polynomial(array(), indeterminate(), galois());
  };

  /**
   * @returns {polynomial} the reciprocal polynomial.
   */
  reciprocal = function () {
    return polynomial(array().reverse(), indeterminate(), galois());
  };

  /**
   * @returns {Number} the degree of the polynomial.
   */
  degree = function () {
    return array().length - 1;
  };

  /**
   * @param {Number} n the index.
   * @returns {Number} the coefficient of the term with exponent `n`.
   */
  coefficient = function (n) {
    return array()[n];
  };

  /**
   * @returns {Number} the leading coefficient of the polynomial.
   */
  leading = function () {
    return coefficient(degree());
  };

  /**
   * @returns {Boolean} the polynomial is monic.
   */
  isMonic = function () {
    return leading() === 1;
  };

  /**
   * @returns {Boolean} the polynomial is *0*.
   */
  isZero = function () {
    return leading() === 0;
  };

  /**
   * @returns {String} a formatted string of the polynomial, e.g., `GF(2): x ^ 2 + x + 1`.
   */
  print = function () {
    var str, i;
    str = "GF(" + galois() + "): ";
    i = degree();
    while (i >= 2) {
      if (coefficient(i) !== 0) {
        str += " + " + coefficient(i) + " " + indeterminate() + " ^ " + i;
      }
      i -= 1;
    }
    if (i === 1 && coefficient(i) !== 0) {
      str += " + " + coefficient(i) + " " + indeterminate();
    }
    str += " + " + coefficient(0);
    str = str.replace(/-/g, "- ");
    str = str.replace(/\+ -/g, "-");
    str = str.replace(/ 1 /g, " ");
    str = str.replace(":  + ", ": ");
    str = str.replace(":  - ", ": - ");
    str = str.replace(" + 0", "");
    return str;
  };

  /**
   * @returns {polynomial} the formal derivative of the polynomial.
   */
  differentiate = function () {
    var outputArray, i;
    outputArray = [];
    for (i = 1; i <= degree(); i += 1) {
      outputArray[i - 1] = modulo(i * coefficient(i), galois());
    }
    return polynomial(outputArray, indeterminate(), galois());
  };

  // Only differentiation; integration is not well-defined

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {Object} the value with key `error` shows if the polynomials have different indeterminates or different orders of fields.
   */
  // Some binary operation can be done only if both operands have the same indeterminate and order of field
  isBinaryOperable = function (inputPolynomial2) {
    var error = undefined;
    if (inputPolynomial2.indeterminate() !== indeterminate()) {
      error = "Unmatched indeterminates";
    }
    if (inputPolynomial2.galois() !== galois()) {
      error = "Unmatched orders of fields";
    }
    return {error:error};
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial
   * @returns {polynomial} the polynomial obtained when `inputPolynomial2` is added to the polynomial, if possible. Otherwise, returns an object with key `error`.
   */
  add = function (inputPolynomial2) {
    var outputArray, i;
    if (typeof isBinaryOperable(inputPolynomial2).error !== "undefined") {
      return isBinaryOperable(inputPolynomial2);
    }
    outputArray = zeros(Math.max(degree(), inputPolynomial2.degree()));
    for (i = 0; i <= degree(); i += 1) {
      outputArray[i] += coefficient(i);
    }
    for (i = 0; i <= inputPolynomial2.degree(); i += 1) {
      outputArray[i] += inputPolynomial2.coefficient(i);
    }
    return polynomial(outputArray, indeterminate(), galois());
  }

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {polynomial} the polynomial obtained when `inputPolynomial2` is subtracted from the polynomial, if possible. Otherwise, returns an object with key `error`.
   */
  subtract = function (inputPolynomial2) {
    var outputArray, i;
    if (typeof isBinaryOperable(inputPolynomial2).error !== "undefined") {
      return isBinaryOperable(inputPolynomial2);
    }
    outputArray = zeros(Math.max(degree(), inputPolynomial2.degree()));
    for (i = 0; i <= degree(); i += 1) {
      outputArray[i] += coefficient(i);
    }
    for (i = 0; i <= inputPolynomial2.degree(); i += 1) {
      outputArray[i] -= inputPolynomial2.coefficient(i);
    }
    return polynomial(outputArray, indeterminate(), galois());
  }

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {polynomial} the polynomial obtained when `inputPolynomial2` is multiplied to the polynomial, if possible. Otherwise, returns an object with key `error`.
   */
  multiply = function (inputPolynomial2) {
    var outputArray, i, j;
    if (typeof isBinaryOperable(inputPolynomial2).error !== "undefined") {
      return isBinaryOperable(inputPolynomial2);
    }
    outputArray = zeros(degree() + inputPolynomial2.degree());
    for (i = 0; i <= degree(); i += 1) {
      for (j = 0; j <= inputPolynomial2.degree(); j += 1) {
        outputArray[i + j] += coefficient(i) * inputPolynomial2.coefficient(j);
      }
    }
    return polynomial(outputArray, indeterminate(), galois());
  }

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {polynomial} the product of `n` repetitions of the polynomial, if possible. Otherwise, returns an object with key `error`.
   */
  power = function (n) {
    var outputPolynomial, i;
    // Power must be non-negative integer
    if (n !== parseInt(n, 10)) {
      return {error: "Power is not an integer"};
    }
    if (n < 0) {
      return {error: "Power is negative"};
    }
    outputPolynomial = polynomial([1], indeterminate(), galois());
    for (i = 1; i <= n ; i += 1) {
      outputPolynomial = outputPolynomial.multiply(clone());
    }
    return outputPolynomial;
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {polynomial} the *largest* polynomial which product with `inputPolynomial2` is *smaller than* or *equal to* the polynomial, if possible. Otherwise, returns an object with key `error`.
   */
  givesQuotient = function (inputPolynomial2) {
    var result;
    result = dividedBy(inputPolynomial2);
    if (typeof result.error === "undefined") {
      return result[0];
    } else {
      return result;
    }
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {polynomial} the *smallest* non-negative polynomial obtained when `inputPolynomial2` and its products with powers of *x* are repeatedly subtracted from the polynomial, if possible. Otherwise, returns an object with key `error`.
   */
  givesRemainder = function (inputPolynomial2) {
    var result;
    result = dividedBy(inputPolynomial2);
    if (typeof result.error === "undefined") {
      return result[1];
    } else {
      return result;
    }
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {polynomial} the polynomial obtained when `inputPolynomial2` is substituted at every instance of the indeterminate. The substituted indeterminate is evaluated in the original field.
   */
  // Composition retains order of finite field but uses new indeterminate
  of = function (inputPolynomial2) {
    var outputPolynomial, i;
    outputPolynomial = polynomial([0], inputPolynomial2.indeterminate(), galois());
    for (i = 0; i <= degree(); i += 1) {
      outputPolynomial = outputPolynomial.add(polynomial(inputPolynomial2.array(), inputPolynomial2.indeterminate(), galois()).power(i).multiply(polynomial(coefficient(i), inputPolynomial2.indeterminate(), galois())));
    }
    return outputPolynomial;
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {Number} the smallest positive integer *k* where the `inputPolynomial2` multiplied by *k* gives *0*, if possible.
   */
  // Get additive order of inputPolynomial2 modulo the current polynomial
  additiveOrderOf = function (inputPolynomial2) {
    var sum, i;
    // Additive order is undefined for infinite fields
    if (!isFinite(galois())) {
      return;
    }
    sum = polynomial(inputPolynomial2.array(), inputPolynomial2.indeterminate(), galois()).givesRemainder(polynomial(array(), inputPolynomial2.indeterminate(), galois()));
    i = 1;
    while(!sum.isZero()) {
      sum = sum.add(polynomial(inputPolynomial2.array(), inputPolynomial2.indeterminate(), galois())).givesRemainder(polynomial(array(), inputPolynomial2.indeterminate(), galois()));
      i += 1;
    }
    return i;
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {Number} the smallest positive integer *k* where `inputPolynomial2` raised to the power of *k* gives *1*, if possible.
   */
  // Get multiplicative order of inputPolynomial2 modulo the current polynomial
  multiplicativeOrderOf = function (inputPolynomial2) {
    var prod, i;
    // Multiplicative order is undefined for infinite fields
    if (!isFinite(galois())) {
      return;
    }
    prod = polynomial(inputPolynomial2.array(), inputPolynomial2.indeterminate(), galois()).givesRemainder(polynomial(array(), inputPolynomial2.indeterminate(), galois()));
    // Multiplicative order is undefined for 0
    if (prod.isZero()) {
      return;
    }
    i = 1;
    while(prod.degree() !== 0 || !prod.isMonic()) {
      prod = prod.multiply(polynomial(inputPolynomial2.array(), inputPolynomial2.indeterminate(), galois())).givesRemainder(polynomial(array(), inputPolynomial2.indeterminate(), galois()));
      i += 1;
      if (i > Math.pow(galois(), degree())) {
        return;
      }
    }
    return i;
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {Boolean} `inputPolynomial2` is a primitive element of the polynomial
   */
  hasPrimitiveElement = function (inputPolynomial2) {
    return multiplicativeOrderOf(inputPolynomial2) === Math.pow(galois(), degree()) - 1;
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {Boolean} `inputPolynomial2` is a factor of the polynomial
   */
  hasFactor = function (inputPolynomial2) {
    return givesRemainder(inputPolynomial2).isZero();
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {Boolean} the polynomial is *0* when evaluated at `inputPolynomial2`.
   */
  hasRoot = function (inputPolynomial2) {
    return of(inputPolynomial2).givesRemainder(polynomial(array(), inputPolynomial2.indeterminate(), galois())).isZero();
  };

  /**
   * @returns {polynomial} the smallest polynomial which is *greater than* the polynomial.
   */
  successor = function () {
    var candidate, i;
    candidate = array();
    candidate[degree() + 1] = 0;
    for (i = 0; i <= degree() + 1; i += 1) {
      candidate[i] = modulo(candidate[i] + 1, galois());
      if (candidate[i] !== 0) {
        return polynomial(candidate, indeterminate(), galois());
      }
    }
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {Boolean} `inputPolynomial2` can be obtained by finding the recursive successor of the polynomial.
   */
  isLessThan = function (inputPolynomial2) {
    var i;
    if (typeof isBinaryOperable(inputPolynomial2).error !== "undefined") {
      return isBinaryOperable(inputPolynomial2);
    }
    if (degree() > inputPolynomial2.degree()) {
      return false;
    }
    if (degree() < inputPolynomial2.degree()) {
      return true;
    }
    for (i = degree(); i >= 0; i -= 1) {
      if (coefficient(i) < inputPolynomial2.coefficient(i)) {
        return true;
      }
      if (coefficient(i) > inputPolynomial2.coefficient(i)) {
        return false;
      }
    }
    return false;
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {Boolean} the polynomial can be obtained by finding the recursive successor of `inputPolynomial2`.
   */
  isGreaterThan = function (inputPolynomial2) {
    var i;
    if (typeof isBinaryOperable(inputPolynomial2).error !== "undefined") {
      return isBinaryOperable(inputPolynomial2);
    }
    if (degree() < inputPolynomial2.degree()) {
      return false;
    }
    if (degree() > inputPolynomial2.degree()) {
      return true;
    }
    for (i = degree(); i >= 0; i -= 1) {
      if (coefficient(i) > inputPolynomial2.coefficient(i)) {
        return true;
      }
      if (coefficient(i) < inputPolynomial2.coefficient(i)) {
        return false;
      }
    }
    return false;
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {Boolean} the polynomial is neither "less than" nor "greater than" `inputPolynomial2`.
   */
  isEqualTo = function (inputPolynomial2) {
    if (typeof isBinaryOperable(inputPolynomial2).error !== "undefined") {
      return isBinaryOperable(inputPolynomial2);
    }
    return !isLessThan(inputPolynomial2) && !isGreaterThan(inputPolynomial2);
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {Array} the additive group of `inputPolynomial2` modulo the polynomial.
   */
  additiveGroupOf = function (inputPolynomial2) {
    var generator, generated, generatedArray;
    generator = polynomial(inputPolynomial2.array(), inputPolynomial2.indeterminate(), galois()).givesRemainder(polynomial(array(), inputPolynomial2.indeterminate(), galois()));
    generated = generator.clone();
    generatedArray = [generated];
    generated = generated.add(generator).givesRemainder(polynomial(array(), inputPolynomial2.indeterminate(), galois()));
    while (!generated.isEqualTo(generator)) {
      generatedArray.push(generated);
      generated = generated.add(generator).givesRemainder(polynomial(array(), inputPolynomial2.indeterminate(), galois()));
      i += 1;
    }
    return generatedArray;
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {Array} the multiplicative group of `inputPolynomial2` modulo the polynomial.
   */
  multiplicativeGroupOf = function (inputPolynomial2) {
    var generator, generated, generatedArray;
    generator = polynomial(inputPolynomial2.array(), inputPolynomial2.indeterminate(), galois()).givesRemainder(polynomial(array(), inputPolynomial2.indeterminate(), galois()));
    generated = generator.clone();
    generatedArray = [generated];
    generated = generated.multiply(generator).givesRemainder(polynomial(array(), inputPolynomial2.indeterminate(), galois()));
    while (!generated.isEqualTo(generator)) {
      generatedArray.push(generated);
      generated = generated.multiply(generator).givesRemainder(polynomial(array(), inputPolynomial2.indeterminate(), galois()));
      i += 1;
    }
    return generatedArray;
  };

  /**
   * @returns {Array} the primitive elements of the polynomial.
   */
  primitiveElements = function () {
    var primitivesArray, primitive;
    if (!isFinite(galois())) {
      return;
    }
    if (isZero()) {
      return;
    }
    primitivesArray = [];
    primitive = polynomial([0, 1], indeterminate(), galois());
    while (degree() > primitive.degree()) {
      if (hasPrimitiveElement(primitive)) {
        primitivesArray.push(primitive.clone());
      }
      primitive = primitive.successor();
    }
    return primitivesArray;
  };

  /**
   * @returns {Array} the factors of the polynomial. The polynoial 0 has infinitely many factors.
   * * For infinite fields, returns
   *   * "undefined" for non-constant polynomials,
   *   * prime factorization for constant polynomials larger than *1*,
   *   * "undefined" for *1* and
   *   * "undefined" for *0*.
   * * For finite fields, returns
   *   * irreducible polynomials factorization for constant polynomials with positive degree,
   *   * "undefined" for constant polynomials larger than *1*,
   *   * empty array for *1* and 
   *   * "undefined" for *0*.
   */
  factors = function () {
    var dividend, factorsArray, factor;
    factorsArray = [];
    if (isZero()) {
      return;
    }
    if (!isFinite(galois())) {
      if (degree() !== 0) {
        return;
      }
      factor = polynomial([2], indeterminate(), galois());
    } else {
      factor = polynomial([0, 1], indeterminate(), galois());
    }
    dividend = clone();
    while (!factor.isGreaterThan(dividend)){
      if (dividend.hasFactor(factor)) {
        factorsArray.push(factor);
        dividend = dividend.givesQuotient(factor);
      } else {
        factor = factor.successor();
      }
    }
    if (isFinite(galois()) && !isMonic()) {
      factorsArray.unshift(polynomial(leading(), indeterminate(), galois()));
    }
    return factorsArray;
  };

  // hasRoots is not defined because there are infinitely many roots in infinitely many polynomials

  /**
   * @returns {Boolean} the polynomial is irreducible (but possibly non-monic).
   */
  isIrreducible = function () {
    return isFinite(galois()) && degree() > 0 && (isMonic() ? factors().length === 1 : factors().length === 1);
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {Boolean} the polynomial in a field is the minimal polynomial of `inputPolynomial2` in the extension field.
   */
  // The polynomial in a field is the minimal polynomial of inputPolynomial2 in the extension field if it is monic irreducible and has inputPolynomial2 as a root
  isMinimalPolynomialOf = function (inputPolynomial2) {
    return isMonic() && isIrreducible() && hasRoot(inputPolynomial2);
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {Array} the primitive polynomials corresponding to the primitive elements returned by `.primitiveElements()`.
   */
  primitivePolynomialsInExtensionOf = function (inputPolynomial2) {
    var candidates, i, primitivePolynomials, primitivePolynomial;
    candidates = primitiveElements();
    primitivePolynomials = [];
    for (i = 0; i < candidates.length; i += 1) {
      primitivePolynomial = polynomial([0, 1], inputPolynomial2.indeterminate(), inputPolynomial2.galois());
      while (!primitivePolynomial.isMinimalPolynomialOf(candidates[i])) {
        primitivePolynomial = primitivePolynomial.successor();
      }
      primitivePolynomials.push(primitivePolynomial);
    }
    return primitivePolynomials;
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {Array} the serialisation counting the occurences of a factor starting with the *0* followed recursively by successors.
   */
  standard = function () {
    var factorsArray, standardArray, index, p, i;
    factorsArray = factors();
    standardArray = [0];
    index = 1;
    p = polynomial([1], indeterminate(), galois());
    i = 0;
    while (i < factorsArray.length) {
      while (p.isLessThan(factorsArray[i])) {
        standardArray[index] = 0;
        p = p.successor();
        index += 1;
      }
      standardArray[index] = 0;
      while (i < factorsArray.length && p.isEqualTo(factorsArray[i])) {
        standardArray[index] += 1;
        i += 1;
      }
      p = p.successor();
      index += 1;
    }
    return standardArray;
  };

  /**
   * @param {Array} standardArray a standard array.
   * @returns {polynomial} the deserialisation of `standardArray`, taking `indeterminate()` and `galois()` from the polynomial.
   */
  fromStandard = function (standardArray) {
    var outputPolynomial, p, i;
    outputPolynomial = polynomial([1], indeterminate(), galois());
    p = polynomial([1], indeterminate(), galois());
    for (i = 2; i < standardArray.length; i += 1) {
      p = p.successor();
      outputPolynomial = outputPolynomial.multiply(p.power(standardArray[i]));
    }
    return outputPolynomial;
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {polynomial} the greatest common divisor of the polynomial and `inputPolynomial2`.
   */
  gcd = function (inputPolynomial2) {
    var standardArray, standardArray2, outputArray, i;
    if (typeof isBinaryOperable(inputPolynomial2).error !== "undefined") {
      return isBinaryOperable(inputPolynomial2);
    }
    standardArray = standard();
    standardArray2 = inputPolynomial2.standard();
    outputArray = [];
    for (i = 0; i < Math.min(standardArray.length, standardArray2.length); i += 1) {
      outputArray[i] = Math.min(standardArray[i], standardArray2[i]);
    }
    return fromStandard(outputArray);
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {polynomial} the lowest common multiple of the polynomial and `inputPolynomial2`.
   */
  lcm = function (inputPolynomial2) {
    var standardArray, standardArray2, outputArray, i;
    if (typeof isBinaryOperable(inputPolynomial2).error !== "undefined") {
      return isBinaryOperable(inputPolynomial2);
    }
    standardArray = standard();
    standardArray2 = inputPolynomial2.standard();
    outputArray = [];
    for (i = 0; i < Math.min(standardArray.length, standardArray2.length); i += 1) {
      outputArray[i] = Math.max(standardArray[i], standardArray2[i]);
    }
    if (typeof standardArray[i] === "undefined") {
      for (i = i; i < standardArray2.length; i += 1) {
        outputArray[i] = standardArray2[i];
      }
    }
    if (typeof standardArray2[i] === "undefined") {
      for (i = i; i < standardArray.length; i += 1) {
        outputArray[i] = standardArray[i];
      }
    }
    return fromStandard(outputArray);
  };

  /**
   * @param {polynomial} inputPolynomial2 another polynomial.
   * @returns {Array} the polynomials *p, xp, xxp, ...*, where *p* is `inputPolynomial2` and *x* is once the indeterminate of the polynomial. The degree of the last element is *1* less than the degree of the polynomial.
   */
  // Used to generate basis of cyclic code
  cycles = function (inputPolynomial2) {
    var basis, basisVector;
    basis = [];
    basisVector = polynomial(inputPolynomial2.array(), inputPolynomial2.indeterminate(), galois());
    while (basisVector.degree() <= degree() - 1) {
      basis.push(basisVector);
      basisVector = basisVector.multiply(polynomial([0, 1], inputPolynomial2.indeterminate(), galois()));
    }
    return basis;
  };

  /**
   * @param {Number} codeLength the length of the desired code.
   * @param {Number} codeDimension the dimension of the desired code.
   * @returns {Array} the arrays corresponding to *a* basis of a cyclic code with code length `codeLength` and code dimension `codeDimension`, taking galois() from the polynomial.
   */
  // Generate polynomials corresponding to a basis of a cyclic code; takes indeterminate() and galois() from the polynomial
  cyclicCodeBasis = function (codeLength, codeDimension) {
    var template, i, unity, minimals, dimension, search, initArray, minimalsChosen, minimal, polynomialsArray, outputArray, j;
    // Code length must be a positive integer
    if (codeLength !== parseInt(codeLength, 10)) {
      throw "Code length is not an integer";
    }
    if (codeLength < 1) {
      throw "Code length is non-positive";
    }
    // Code dimension must be a non-negative integer smaller than code length
    if (codeDimension !== parseInt(codeDimension, 10)) {
      throw "Code dimension is not an integer";
    }
    if (codeDimension < 0) {
      throw "Code dimension is negative";
    }
    if (codeDimension >= codeLength) {
      throw "Code dimension is not smaller than code length";
    }
    template = [-1];
    for (i = 1; i < codeLength ; i += 1) {
      template.push(0);
    }
    template.push(1);
    unity = polynomial(template, indeterminate(), galois());
    minimals = unity.factors();
    if (typeof minimals === "undefined") {
      return;
    }
    dimension = [];
    for (i = 0; i < minimals.length; i += 1) {
      dimension[i] = minimals[i].degree();
    }
    // Return "undefined" if a cyclic code of length codeLength and dimension codeDimension cannot be found
    search = function (isChosen, sum) {
      if (isChosen.length === dimension.length) {
        return;
      }
      if (sum + dimension[isChosen.length] === codeLength - codeDimension) {
        isChosen.push(true);
        return isChosen;
      }
      if (sum + dimension[isChosen.length] > codeLength - codeDimension) {
        while (!isChosen[isChosen.length - 1]) {
          isChosen.pop();
          if (isChosen.length === 0) {
            return;
          }
        }
        isChosen[isChosen.length - 1] = false;
        sum -= dimension[isChosen.length - 1];
        return search(isChosen, sum);
      }
      sum += dimension[isChosen.length];
      isChosen.push(true);
      return search(isChosen, sum);
    };
    initArray = [];
    minimalsChosen = search(initArray.slice(0), 0);
    while (typeof minimalsChosen === "undefined" && initArray.length <= dimension.length - 2) {
      initArray.push(false);
      minimalsChosen = search(initArray.slice(0), 0);
    }
    if (typeof minimalsChosen === "undefined") {
      return;
    }
    minimal = polynomial([1], indeterminate(), galois());
    for (i = 0; i < minimals.length; i += 1) {
      if (minimalsChosen[i]) {
        minimal = minimal.multiply(minimals[i]);
      }
    }
    polynomialsArray = unity.cycles(minimal);
    outputArray = [];
    for (i = 0; i < polynomialsArray.length; i += 1) {
      outputArray[i] = polynomialsArray[i].array();
      for (j = outputArray[i].length; j < codeLength; j += 1) {
        outputArray[i][j] = 0;
      }
    }
    return outputArray;
  };

  return {
    array: array,
    indeterminate: indeterminate,
    galois: galois,
    clone: clone,
    reciprocal: reciprocal,
    degree: degree,
    coefficient: coefficient,
    leading: leading,
    isMonic: isMonic,
    isZero: isZero,
    print: print,
    differentiate: differentiate,
    isBinaryOperable: isBinaryOperable,
    add: add,
    subtract: subtract,
    multiply: multiply,
    power: power,
    givesQuotient: givesQuotient,
    givesRemainder: givesRemainder,
    of: of,
    additiveOrderOf: additiveOrderOf,
    multiplicativeOrderOf: multiplicativeOrderOf,
    hasPrimitiveElement: hasPrimitiveElement,
    hasFactor: hasFactor,
    hasRoot: hasRoot,
    successor: successor,
    isLessThan: isLessThan,
    isGreaterThan: isGreaterThan,
    isEqualTo: isEqualTo,
    additiveGroupOf: additiveGroupOf,
    multiplicativeGroupOf: multiplicativeGroupOf,
    primitiveElements: primitiveElements,
    factors: factors,
    isIrreducible: isIrreducible,
    isMinimalPolynomialOf: isMinimalPolynomialOf,
    primitivePolynomialsInExtensionOf: primitivePolynomialsInExtensionOf,
    standard: standard,
    fromStandard: fromStandard,
    gcd: gcd,
    lcm: lcm,
    cycles: cycles,
    cyclicCodeBasis: cyclicCodeBasis
  };

} catch (error) {
  return {error: error};
}

};
