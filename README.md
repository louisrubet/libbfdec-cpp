# libbfdec-cpp

C++ class for **Fabrice Bellard libbf library** (decimal part only)

This ["Tiny Big Float library" libbf](https://bellard.org/libbf/) allows to do IEEE 754 floating point calculations (both binary and decimal) with arbitrary precision and transcendent functions. [This online calculator](http://numcalc.com/) uses libbf.

- This class is mostly done to replace [mpreal class](https://github.com/advanpix/mpreal/), which is a C++ helper class over [mpfr](http://mpfr.org),
- for being used in [rpn functional language](https://github.com/louisrubet/rpn),
- libbf is included without modifications as a submodule from [this github repository](https://github.com/rurban/libbf).

Fake `quickjs.h` and `quickjs-config.h` are provided so libbfdec-cpp can be build without including [quickjs](https://bellard.org/quickjs/).

## build

Although [libbf](https://github.com/rurban/libbf) is available for several operating systems, the build was only tested on linux.

Are needed `cmake` and `gcc`. No other development lib is needed.

```
mkdir build
cmake -B build && make -C build
```

## run test

```
build/libbfdec-cpp-test
```
