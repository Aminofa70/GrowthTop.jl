"""
    get_material_matrix(E, ν)

Compute the 4th-order elasticity (stiffness) tensor for a 2D isotropic linear elastic material,
given the Young's modulus `E` and Poisson's ratio `ν`.

# Arguments
- `E`: Young's modulus — measures material stiffness.
- `ν`: Poisson's ratio — the ratio of transverse strain to axial strain.

# Returns
- A `SymmetricTensor{4,2}` representing the full 4th-order stiffness tensor.

# Method
1. Constructs the 2D stiffness matrix in Voigt notation:
   C_voigt = (E / (1 - ν^2)) * [1.0  ν     0.0
                                ν    1.0   0.0
                                0.0  0.0  (1 - ν)/2]

2. Converts the Voigt matrix to a full 4th-order symmetric tensor.
"""
function get_material_matrix(E, ν)
    C_voigt = E * [1.0 ν 0.0; ν 1.0 0.0; 0.0 0.0 (1 - ν)/2] / (1 - ν^2)
    return fromvoigt(SymmetricTensor{4,2}, C_voigt)
end






