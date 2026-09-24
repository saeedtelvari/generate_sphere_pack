# Post-Depositional Diagenesis Engine (Compaction & Cementation)

Type: grilling
Status: resolved
Blocked by: 01

## Question

How should mechanical compaction (anisotropic vertical strain tensor) and cementation (Euclidean distance transform dilation vs curvature-driven overgrowth) be parameterized to mirror realistic sandstone diagenesis?

Specifically:
1. In mechanical compaction, should the vertical strain $\varepsilon_{zz}$ be applied directly to grain centers (introducing interpenetration and contact squashing), or should grain deformation / flattening be modeled geometrically?
2. In cementation, how should secondary mineral precipitation (quartz overgrowths, pore-filling calcite) be modeled? Does a thresholded Euclidean Distance Transform (EDT) on the grain phase realistically reproduce preferential filling of small pore throats before larger pore bodies?
3. How should secondary dissolution (leaching of grains or cements to create secondary micro-porosity) be exposed in the API?

## Answer

### 1. Mechanical Compaction (Lagrangian Uniaxial Strain)
- **Transformation**: Compaction is executed directly on `GrainPack` via vertical uniaxial strain $\varepsilon_{zz} \in [0, 0.25]$:
  $$z' = (1 - \varepsilon_{zz}) z, \quad L_z' = (1 - \varepsilon_{zz}) L_z$$
- **Physical Effect**: Squeezes grain centroids together, creating realistic grain interpenetration and contact indentation observed in thin sections of deeply buried sandstones. Porosity drops from the uncompacted deposition baseline ($\sim 0.40$) down to $\sim 0.28 - 0.34$, while establishing vertical-horizontal permeability anisotropy ($k_v < k_h$).

### 2. Cementation / Overgrowth (Eulerian Distance Transform Dilation)
- **Mechanism**: Evaluates the 3D Euclidean Distance Transform (EDT) into the pore phase:
  $$D(\mathbf{x}) = \text{EDT}(\text{matrix} == \text{Pore})$$
- **Throat Filling Dynamics**: In a porous rock, narrow pore throats have small clearance to surrounding solid surfaces compared to large pore bodies. Applying a cementation threshold $d_{cement}$ (`matrix[D <= d_cement] = True`):
  - Completely seals and occludes tight pore throats.
  - Leaves large pore bodies partially open with coated grain rims.
  - Accurately captures syntaxial quartz overgrowths and poikilotopic calcite cementation, reducing porosity from $0.30$ down to $0.05 - 0.20$.

### 3. Secondary Porosity & Leaching
- Secondary dissolution is parameterized via `dissolve(fraction, mode='matrix'|'grain_cores')`:
  - `matrix`: Random cluster erosion expanding pore boundaries.
  - `grain_cores`: Hollows out internal grain cores to simulate unstable feldspar leaching in arkosic sandstones.
