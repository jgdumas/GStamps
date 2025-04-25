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
           `sudo apt install git make g++ libboost-dev`
           (and optionally `sudo apt install pkg-config libgmp-dev libgivaro-dev`).
- Then just run `make`, in order to produce the following executable programs
- See also [`bin/auto-docker.run`](https://github.com/jgdumas/gstamps/blob/main/bin/auto-docker.run)


**Programs**:
|  |  |
| :--------- | :------ |
|`bin/cover`| LPSP: Computes the cover of a basis with s stamps |
|`bin/basis`| GPSP: Divide & Conquer basis computation of k denominations and s stamps |
|`bin/dynprg`| GPSP: Dynamic programming search of best Divide & Conquer cut |
|  |  |


**Search Tools**:
|  |  |
| :--------- | :------ |
|`bin/search`| Smallest basis reaching N with s stamps |
|`bin/complement`| Exhausts all additional denominations in parallel |
|`bin/supplement`| Exhausts additional denominations several values at a time |
|`bin/brute`| Exhaustive search of an extremal basis |
|  |  |

**Other basis**:
|  |  |
| :--------- | :------ |
|`bin/fibo`| Fibonacci basis |
|`bin/geom`| Geometric progression basis |
|`bin/ab`| Alter & Barnett basis |
|  |  |

**Usage**:
- #k: k denominations in the basis
- #s: the basis is for s stamps
- #v: verbosisty level
- #r: after r recursive levels (rlevel), stops searching for the best cut, just use the midpoint
- #a: if true provides only a lower bound on the cover of the basis (approximate), otherwise computes the cover exactly

**Examples**:
- `./bin/basis 4 2`: produces a basis of 4 denominations for 2 stamps (1 3 5 6, attaining all integers 1..12)
- `echo '1 3 5 6' | ./bin/cover 2`: the basis can cover (all integers up to 12) with 2 stamps
- `./bin/search 4 1024 1`: produces a basis (with 14 denominations), that can cover at least all integers up to 1024 with 4 stamps