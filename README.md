--------------------------------------------------------------------------------
# GStamps: C++ routines for the Global Postage Stamp Problem
--------------------------------------------------------------------------------

The postage stamp problem states that an envelope may be franked with
a total of at most s stamps while one has available an (additive)
basis of k stamp denominations (integers: a1 < a2 < ... < ak).

- The local postage stamp problem (LPSP) is the determination of the
smallest integer not attainable by a given basis with at most s stamps
- The global postage stamp problem (GPSP) is the determination of a
basis with the largest LPSP for given parameters k and s.

**Authors**: 
Jean-Guillaume Dumas,
Alexis Galan,
Bruno Grenet,
Aude Maignan,
Daniel S. Roche


**Requirements**:
- C++
- Optional: [Givaro](https://github.com/linbox-team/givaro), dev: headers & library (version ≥ 4.2.0)


**Installation**:
- Requires some distribution packages like:
           `sudo apt install make g++`
           (and optionally `sudo apt install  pkg-config libgivaro-dev).
- Then just run `make`, in order to produce the following executable programs
- See also [`bin/auto-docker.run`](https://github.com/jgdumas/gstamps/blob/main/bin/auto-docker.run)


**Programs**:
|  |  |
| :--------- | :------ |
|`bin/cover`| LPSP: Computes the cover of a basis with s stamps |
|`bin/basis`| GPSP: Direct basis computation of k denominations and s stamps |
|`bin/dynprg`| GPSP: Dynamic programming search of best Divide & Conquer cut |
|`bin/search`| Smallest basis reaching N with s stamps |
|`bin/supplement`  | Exhausts additional denominations one value at a time |
|`bin/complement`| Exhausts all additional denominations in parallel |
|`bin/brute`| Exhaustive search of an extremal basis |
|  |  |


**Tools**:
|  |  |
| :--------- | :------ |
|`bin/fibo`| Fibonacci basis |
|`bin/geom`| Geometric progression basis |
|`bin/ab`| Alter & Barnett basis |
|  |  |


**Examples**:
- `./bin/basis 4 2`: produces the basis "1 3 5 6" attaining all integers 1..12
- `./bin/`: the 
