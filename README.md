# libbfdec-cpp

C++ class for **Fabrice Bellard libbf library** (decimal part only)

- this work is mostly done to replace [mpreal class](https://github.com/advanpix/mpreal/) (C++ helper class over [mpfr](http://mpfr.org))
- for being used in [rpn functional language](https://github.com/louisrubet/rpn)
- Fabrice Bellard's libbf is taken without modifications in [this github repository](https://github.com/rurban/libbf).

Fake `quickjs.h` and `quickjs-config.h` are provided so libbfdec-cpp can be build without including [quickjs](https://bellard.org/quickjs/).

## build

Although [libbf](https://github.com/rurban/libbf) is available for several operating systems, the build was only tested on linux.

Are needed `cmake` and `gcc`. No other development lib is needed.

```
mkdir build
cmake -B build && make -C build
```

## run est

```
build/libbfdec-cpp-test
```
