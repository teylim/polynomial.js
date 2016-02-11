# polynomial.js
***a JavaScript library for working with polynomials over finite fields***

License is omitted for the time being, until someone finds the files useful.

Documentation is found in the **output/** folder. It is generated using [jsdox.js](http://jsdox.org/) licensed under the MIT license.

**polynomial.js** is the main file. Polynomials are created using function declaration `polynomial(inputArray, inputIndeterminate, inputGalois)`, not via `new`.

**tests.html** is used to check for errors using a console.

Important Notes:
* The use of the term "infinite fields" for *{0, 1, 2, ...}* is just a matter of programming convenience. It is mathematically *wrong* because there does not exist an additive inverse or a multiplicative inverse for some elements.
* Define the increase in degree across polynomials as increasing and between polynomials having the same degree, the ordering *0, 1, ...* for each coefficient as increasing.

Copyright 2016 @litena