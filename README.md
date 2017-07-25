ipws - Classify Weight Systems with Interior Point
==================================================

Building
--------

Requirements:

- CMake
- C++14 compiler
- boost library

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
dnf install glibc-static libstdc++-static boost-static
```

Then run:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DLINK_STATICALLY=ON -DDIMENSION=4 -DR_NUMERATOR=1 -DR_DENOMINATOR=2 ..
make
```

Usage
-----

### Print all IP weight system candidates

This does all in one run, and is feasible for `d - 2r < 4`.

```
./ipws --find-candidates --print-ws
```

### Print all IP weight systems

This does all in one run, and is feasible for `d - 2r < 4`.

```
./ipws --find-ip --print-ws
```

### Find all IP weight systems, parallelizable

First, generate a list of weight system pairs and the list of weight systems found along the way:

```
./ipws --find-pairs --pairs-out pairs --ws-out ws
```

Next iterate over all the weight system pairs and generate weight systems:

```
./ipws --find-candidates --start 0 --count 100 --pairs-in pairs --ws-out ws-0
./ipws --find-candidates --start 100 --count 100 --pairs-in pairs --ws-out ws-100
...
```

Combine all the so found weight systems into one big list, two at a time:

```
./ipws --combine-ws ws-0 --ws-in ws-100 --ws-out ws-a
...
```
