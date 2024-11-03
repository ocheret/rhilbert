# rhilbert

Package rhilbert provides functions to convert back and forth between
Hilbert Curve values and points in n-space.

This uses an extension to the algorithm described in:

```
"Alternative Algorithm for Hillbert's Space-Filling Curve," A R
Butz, IEEE Trans. on Computers, April 1971, pp 424-426.
```
*(NOTE: Butz was a rabid Holocaust denier but this does not mean we can't use this work.)*

The algorithm extension added here allows for conversion back and forth and not just in one direction.

For dimensionality n and order m:
- The Hilbert number varies from 0 to H = (2 ** (n * m)) - 1
- The integer grid varies from 0 to G = (2 ** m) - 1 along each axis

Examples:
- n = 2, m = 4 -> H = 255, G = 15
- n = 3, m = 5 -> H = 32767, G = 31

For n = 2, m = 2 there are H=15 unique Hilbert numbers and the X
and Y coordinates vary between 0 and G=3.

The correspondence between Hilbert numbers and coordinates is shown
here for n = 2 and m = 2:
```text
    3| 5---6   9--10
     | |   |   |   |
    2| 4   7---8  11
   Y | |           |
    1| 3---2  13--12
     |     |   |
    0| 0---1  14--15
     +--------------
       0   1   2   3
           X
```
```text
   Hilbert    Coordinate
   Number      [X Y]
   -------    ----------
     0   <->   [0 0]
     1   <->   [1 0]
     2   <->   [1 1]
     3   <->   [0 1]
     4   <->   [0 2]
     5   <->   [0 3]
     6   <->   [1 3]
     7   <->   [1 2]
     8   <->   [2 2]
     9   <->   [2 3]
    10   <->   [3 3]
    11   <->   [3 2]
    12   <->   [3 1]
    13   <->   [2 1]
    14   <->   [2 0]
    15   <->   [3 0]
```
