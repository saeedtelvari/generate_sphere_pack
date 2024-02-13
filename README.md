# Generate 3D randomly generated sphere packs (porous media)

Features
- __Brute-force packing__: Implements a simple algorithm for randomly placing spheres while checking for collisions/overlapping.
- __Configurable parameters__: Able to control the number of spheres, size distribution, and bounding box size.
- __Output formats__: Export the generated pack as a csv file containing sphere centers (x, y, z) and radii (r), which can be converted into a 3D image using porespy.
------------------
**Getting Started:**
1. **Clone the repository:**
```
git clone https://github.com/saeedtelvari/generate_sphere_pack.git
```
2. **Install dependencies:**
```
pip install -r requirements.txt
```
3. **Import class and use its functions:**
```
from genrandsp import *
```
------------------
**Notes:**
* An example notebook demonstrating the usage of this repository is available.
* The functional programming approach initially implemented in function.py was later transformed into object-oriented programming (OOP) as it proved to be a more suitable paradigm.
------------------
**Contributions:**
* Feel free to reach out if you're interested in contributing by fixing bugs or implementing new features and algorithms. We welcome your contributions!
* Brace yourself for some wild variable naming and messy code ahead! :))))
