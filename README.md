ipws - Classify Weight Systems with Interior Point
==================================================

Building
--------

Requirements:

- CMake >= 3.1
- C++14 compiler
- boost

The dimension `d` and the index `r` of the weight systems are determined at compile time.
For an ordinary release build with `d=4` and `r=1/2` run:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DDIMENSION=4 -DR_NUMERATOR=1 -DR_DENOMINATOR=2 ..
make
```

### Link Statically

If the standard library or boost is too old on the machine where the binary is supposed to run, it can be statically linked instead.
Install the required packages (on Fedora):

```
dnf install glibc-static libstdc++-static
```

Then run:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DLINK_STATICALLY=ON -DDIMENSION=4 -DR_NUMERATOR=1 -DR_DENOMINATOR=2 ..
make
```
