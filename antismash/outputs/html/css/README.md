antiSMASH styling instructions
==============================

Please use the `lessjs` files for styling, don't edit the `css` files directly.

Requirements
------------

This requires the `npm` utility installed, and the `less` and `less-plugin-clean-css` packages.
Once you have `npm` installed, you can install the latter by running:

```
npm install -g less less-plugin-clean-css
```

Generating new `css` files
--------------------------

Provided you have `lessc` and `less-plugin-clean-css` installed in your path, simply run `make`
in the directory containing the `.less` files and this `README.md`.
