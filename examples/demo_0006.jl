using Revise
using Ferrite
using PUPM
using GrowthTop
# done!
### topology optmization three points bending Tone
# Function to create the grid
par = DynamicParams()
function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly))
    ]
    grid = Ferrite.generate_grid(Ferrite.Quadrilateral, (nx, ny), corners)
    addnodeset!(grid, "clamped", x -> x[1] ≈ 0.0  && x[2] ≈ 0.0) #fixed in x-ydirection
    addnodeset!(grid, "support_2", x -> x[1] ≈ Lx && x[2] ≈ 0.0) # fixed in y direction
    addnodeset!(grid, "nodal_force", x -> x[1] ≈ Lx/2 && x[2] ≈ Ly) # nodal_force
    return grid
end

# Function to create CellValues and FacetValues
function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{Ferrite.RefQuadrilateral, order}()^dim
    qr = Ferrite.QuadratureRule{Ferrite.RefQuadrilateral}(2)
    qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefQuadrilateral}(1)
    cell_values = Ferrite.CellValues(qr, ip)
    facet_values = Ferrite.FacetValues(qr_face, ip)
    return cell_values, facet_values
end

# Function to create DofHandler
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()^2)
    Ferrite.close!(dh)
    return dh
end

# Function to create Dirichlet boundary conditions
function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getnodeset(dh.grid, "clamped"), (x, t) -> [0.0, 0.0], [1, 2]))
    add!(ch, Dirichlet(:u, getnodeset(dh.grid, "support_2"), (x, t) -> [0.0], [2]))
    Ferrite.close!(ch)
    return ch
end

# Define parameters for the plate and mesh
Lx, Ly = 2.0, 1.0  # Plate dimensions
nx, ny = 120, 60   # Number of elements along x and y

grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid
par.tnele = length(grid.cells)  # Total number of elements

par.grid = grid
# Create DOF handler and constraints
par.dh = create_dofhandler(grid)
par.ch = create_bc(par.dh)

# Create CellValues and FacetValues
par.cell_values, par.facet_values = create_values()

# Define loads
par.loads = [LoadCondition("nodal_load", [0.0, -1.0])]  # Load applied to the "traction" facet

# Material properties
par.E0 = 1.0                # Initial Young's modulus
E = fill(par.E0, Ferrite.getncells(grid))  # Initialize Young's modulus for all cells

par.ν = 0.3  # Poisson's ratio

# Optimization parameters
par.Emin = 1e-4             # Minimum Young's modulus
par.Emax = 1.0              # Maximum Young's modulus
par.ρ0 = 1.0                # Initial density
par.tol = 1e-2             # Convergence tolerance
# gamma_values = [1, 2, 3]                     # Penalty factor
# eta_values = [π/(2.5) , π /(3.5) , π/(4.0)]               # Filter parameter
# k_values = [1 , 4, 8]                     # Sensitivity parameter
# volfrac_values =[0.0, 0.25, 0.5, 0.75]                   # Volume fraction


γ= 1
η = π/(3.5)
k = 4
volfrac = 0.25
# Neumann BC facet set
par.Neumann_bc = Ferrite.getnodeset(grid, "nodal_force")  # Nodes on the edge

file_name = "optim"
dir = "/Users/aminalibakhshi/Desktop/vtu_file/"
remove_vtk_files(dir)
par.max_itr = 200
top_2d(Tone , par, E, k, γ, η ,volfrac, file_name, dir)