### define the type of topology optimization problem
struct UPM end
struct SIMP end
struct Tone end
struct Ttwo end
################################################################################

function transfer_to_density(Enew::Array{Float64,1}, E0::Float64, ρ0::Float64, γ::Int64)
    # Define the minimum and maximum allowed density values
    ρmin, ρmax = 0.0, 1.0

    # Calculate the density based on the given formula
    ρ = ρ0 .* ((Enew / E0) .^ (1 / γ))

    # Clamp the density values to ensure they stay within the range [ρmin, ρmax]
    ρ = [clamp(x, ρmin, ρmax) for x in ρ]

    # Return the resulting density array
    return ρ
end
################################################################################
function transfer_to_young(ρnew::Array{Float64,1}, E0::Float64,
    ρ0::Float64, γ::Int64, Emin::Float64, Emax::Float64)

    # Calculate Young's modulus using the formula:
    #   E = E0 * ((ρ / ρ0) ^ γ)
    # Element-wise operations are performed on the input array `ρnew`.
    Enew = E0 * (ρnew / ρ0) .^ γ

    # Clamp the Young's modulus values to ensure they stay within the range [Emin, Emax].
    # This ensures no values fall below the minimum or exceed the maximum.
    Enew = [clamp(x, Emin, Emax) for x in Enew]

    # Return the resulting array of clamped Young's modulus values.
    return Enew
end

################################################################################
function filter_density_to_vf!(density, vf, tnele, eta)
    # Define the minimum and maximum density values allowed.
    rhomin, rhomax = 0.01, 1.0

    # Inner function to transform density values based on the current threshold (rhotr).
    function transform(rholoc, rhotr, eta, rhomin, rhomax)
        if rholoc < rhotr
            # Below the threshold, assign minimum density value.
            rhotrans = rhomin
        elseif rholoc > rhotr + 1.0 / tan(eta)
            # Above the threshold + filter range, assign maximum density value.
            rhotrans = rhomax
        else
            # Linearly transform values within the threshold range.
            rhotrans = tan(eta) * (rholoc - rhotr)
        end
        return rhotrans
    end

    # Define the bounds for rhotr (threshold value) based on the filter range.
    rhomaxbound = -1.0 / tan(eta)  # Minimum value for a volume fraction of 0.
    rhominbound = 1.0             # Maximum value for a volume fraction of 1.

    # Initialize the error as a large number for the iteration loop.
    error = 10.0
    rhotr = 0.0  # Initialize the threshold value.

    # Iteratively adjust `rhotr` to achieve the target volume fraction (`vf`).
    while error > 0.001
        rhotr = (rhominbound + rhomaxbound) / 2  # Midpoint of the bounds.

        # Calculate the volume fraction for the lower, upper, and current threshold.
        sumdmin = sum(transform(density[i], rhominbound, eta, rhomin, rhomax) for i in eachindex(density)) / tnele
        sumdmax = sum(transform(density[i], rhomaxbound, eta, rhomin, rhomax) for i in eachindex(density)) / tnele
        sumdmid = sum(transform(density[i], rhotr, eta, rhomin, rhomax) for i in eachindex(density)) / tnele

        # Update bounds based on the relation between the calculated and target volume fractions.
        if (sumdmin - vf) / (sumdmid - vf) > 0
            rhominbound = rhotr
        elseif (sumdmax - vf) / (sumdmid - vf) > 0
            rhomaxbound = rhotr
        else
            println("Problem out of bounds:", sumdmax, sumdmin, sumdmid, vf)
        end

        # Update the error value to track convergence.
        error = abs(vf - sumdmid)
    end

    # Apply the final transformation to the density values using the calculated threshold.
    for i in eachindex(density)
        densloc = transform(density[i], rhotr, eta, rhomin, rhomax)
        density[i] = densloc
    end

    # Return the adjusted density array.
    return density
end
################################################################################
################################################################################
################################################################################
function update_young_UPM(k, E, H, Emax, Emin, E0, γ)
    # Initialize updated Young's modulus array
    Enew = similar(E)

    # Calculate mean and standard deviation of H
    H_mean = mean(H)
    H_std = std(H)

    # Initialize α array
    α = similar(E)

    for i in eachindex(E)
        # Compute part1 as a scaled function of H[i] and E[i]
        part1 = H[i] * (E[i] / E0)^((γ - 1) / γ)

        # Compute normalized adjustment value α[i]
        if part1 - H_mean < 0
            α[i] = -((abs(part1 - H_mean)) / (k * H_std))^γ
        else
            α[i] = ((part1 - H_mean) / (k * H_std))^γ
        end

        # Update E[i] with α adjustment
        α_scalar = α[i]
        Enew[i] = E[i] * (1 + α_scalar)

        # Apply limits to ensure Enew[i] stays within [Emin, Emax]
        Enew[i] = clamp(Enew[i], Emin, Emax)
    end

    return Enew
end
################################################################################
function update_young_SIMP(E, H, Emax, Emin, E0, γ)
    

    # Initialize arrays for updated values and intermediate results
    Enew = similar(E)
    α = similar(E)

    # Iterate over each element to compute updated values
    for i in eachindex(E)
        # Compute the intermediate value α for the current element
        α[i] = (H[i] * (E[i] / E0)^((γ - 1) / γ))^γ

        # Update Young's modulus based on the calculated α value
        Enew[i] = E[i] * (1 + α[i])

        # Enforce the upper and lower bounds on the updated modulus
        Enew[i] = clamp(Enew[i], Emin, Emax)
    end

    return Enew
end
################################################################################
function update_young_Tone(k, E, H, Emax, Emin, E0, γ)
    

    # Initialize updated Young's modulus array
    Enew = similar(E)

    # Calculate mean and standard deviation of H
    H_mean = mean(H)
   # H_std = std(H)

    # Initialize α array
    α = similar(E)

    for i in eachindex(E)
        # Compute part1 as a scaled function of H[i] and E[i]
        part1 = H[i] * (E[i] / E0)^((γ - 1) / γ)

        # Compute normalized adjustment value α[i]
        if part1 - H_mean < 0
            α[i] = -((abs(part1 - H_mean)) / (k ))^γ
        else
            α[i] = ((part1 - H_mean) / (k ))^γ
        end

        # Update E[i] with α adjustment
        α_scalar = α[i]
        Enew[i] = E[i] * (1 + α_scalar)

        # Apply limits to ensure Enew[i] stays within [Emin, Emax]
        Enew[i] = clamp(Enew[i], Emin, Emax)
    end

    return Enew
end
################################################################################
function update_young_Ttwo(E, H, Emax, Emin, E0, γ, B_Δt)
    Enew = similar(E)
    #H_mean = mean(H)

    H_mean = 0.0
    ρ0 = 1.0
    residual = (E0*B_Δt*γ)/(ρ0^2)

    α = similar(E)
    
    for i in eachindex(E)
        # Compute part1 as a scaled function of H[i] and E[i]
        part1 = H[i] * (E[i] / E0)^((γ - 1) / γ)

        # Compute normalized adjustment value α[i]
        if part1 - H_mean < 0
            α[i] = -((abs(part1 - H_mean)))^γ
        else
            α[i] = ((part1 - H_mean))^γ
        end

        # Update E[i] with α adjustment
        α_scalar = α[i]*(residual^γ)
        Enew[i] = E[i] + E0 * α_scalar

        # Apply limits to ensure Enew[i] stays within [Emin, Emax]
        Enew[i] = clamp(Enew[i], Emin, Emax)
    end
   
    return Enew
    
end
################################################################################
function top_2d(::Type{UPM} , par::DynamicParams, E, k, γ, η ,volfrac, name_of_file::String, directory::String)
    grid = par.grid
    dh = par.dh
    #E = par.E
    # nx = par.nx ; ny = par.ny ; nz = par.nz
    tnele = par.tnele
    E0 = par.E0 ; Emin = par.Emin ; Emax = par.Emax
    #k = par.k ; γ = par.γ ; volfrac = par.vf; 
    #η = par.η; 
    ρ0 = par.ρ0
    max_itr = par.max_itr ; tol = par.tol

    loop = 1

    ## Initial FEM solve
    fem = fem_solver_combine(par, E)
    compliance = fem.compliance
    H = fem.H
    W_tot = sum(fem.U)
    strain_energy_vector = [W_tot, W_tot * 10]
    A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]
    println("Iter $loop: C = $compliance, ΔR = $A")
    ## Iterative optimization loop
    while abs(A) > tol && loop <= max_itr
        fem = fem_solver_combine(par, E)
        compliance = fem.compliance
        H = fem.H
        W_tot = sum(fem.U)
        
        if volfrac == 0.0
            Enew = update_young_UPM(k, E, H, Emax, Emin, E0, γ)
            Enew_frac = Enew
        elseif volfrac > 0.0
            Enew = update_young_UPM(k, E, H, Emax, Emin, E0, γ)
            ρ = transfer_to_density(Enew, E0, ρ0, γ)
            ρnew = filter_density_to_vf!(ρ, volfrac, tnele, η)
            Enew_frac = transfer_to_young(ρnew, E0, ρ0, γ, Emin, Emax)
        else
            error("Invalid value for volfrac")
        end

        # Update E in par so that fem_solver uses the updated material distribution
        #par.E = Enew_frac
        E = Enew_frac

        fem = fem_solver_combine(par, E)
        compliance = fem.compliance
        H = fem.H
        W_tot = sum(fem.U)

        strain_energy_vector[1] = strain_energy_vector[2]
        strain_energy_vector[2] = W_tot
        A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]

        loop += 1
        println("Iter $loop: C = $compliance, ΔR = $A")


    end

    ## Handle termination
    if loop > max_itr
        compliance = -1
    end

    # Final write
    fem = fem_solver_combine(par, E)
    compliance = fem.compliance
    u = fem.u
    σ = fem.σ
    ε = fem.ε
    E_node = fem.E_node
    U = fem.U
    full_path = joinpath(directory, name_of_file)

    VTKGridFile(full_path, dh) do vtk
        write_solution(vtk, dh, u)
        for (j, key) in enumerate(("11", "22", "12"))
            write_cell_data(vtk, σ[j], "stress_" * key)
        end
        for (j, key) in enumerate(("11", "22", "12"))
            write_cell_data(vtk, ε[j], "strain_" * key)
        end
        write_cell_data(vtk, E, "Young's modulus")
        write_cell_data(vtk, U, "Strain Energy")
        write_node_data(vtk, E_node, "Nodal Young's modulus")
        Ferrite.write_cellset(vtk, grid)
    end

    return compliance

end
################################################################################
function top_2d(::Type{SIMP} , par::DynamicParams, E, γ, η ,volfrac, name_of_file::String, directory::String)
    grid = par.grid
    dh = par.dh
    #E = par.E
    # nx = par.nx ; ny = par.ny ; nz = par.nz
    tnele = par.tnele
    E0 = par.E0 ; Emin = par.Emin ; Emax = par.Emax
    #γ = par.γ ; volfrac = par.vf; 
    #η = par.η; 
    ρ0 = par.ρ0
    max_itr = par.max_itr ; tol = par.tol

    loop = 1

    ## Initial FEM solve
    fem = fem_solver_combine(par, E)
    compliance = fem.compliance
    H = fem.H
    W_tot = sum(fem.U)
    strain_energy_vector = [W_tot, W_tot * 10]
    A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]
    println("Iter $loop: C = $compliance, ΔR = $A")
    ## Iterative optimization loop
    while abs(A) > tol && loop <= max_itr
        fem = fem_solver_combine(par, E)
        compliance = fem.compliance
        H = fem.H
        W_tot = sum(fem.U)
        # Enew = update_young_SIMP(E, H, Emax, Emin, E0, γ)
        # ρ = transfer_to_density(Enew, E0, ρ0, γ)
        # ρnew = filter_density_to_vf!(ρ, volfrac, tnele, η)
        # Enew_frac = transfer_to_young(ρnew, E0, ρ0, γ, Emin, Emax)
        if volfrac == 0.0
            Enew = update_young_SIMP(E, H, Emax, Emin, E0, γ)
            Enew_frac = Enew
        elseif volfrac > 0.0
            Enew = update_young_SIMP(E, H, Emax, Emin, E0, γ)
            ρ = transfer_to_density(Enew, E0, ρ0, γ)
            ρnew = filter_density_to_vf!(ρ, volfrac, tnele, η)
            Enew_frac = transfer_to_young(ρnew, E0, ρ0, γ, Emin, Emax)
        else
            error("Invalid value for volfrac")
        end

        # Update E in par so that fem_solver uses the updated material distribution
        #par.E = Enew_frac
        E = Enew_frac

        fem = fem_solver_combine(par, E)
        compliance = fem.compliance
        H = fem.H
        W_tot = sum(fem.U)

        strain_energy_vector[1] = strain_energy_vector[2]
        strain_energy_vector[2] = W_tot
        A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]

        loop += 1
        println("Iter $loop: C = $compliance, ΔR = $A")

    end

    ## Handle termination
    if loop > max_itr
        compliance = -1
    end

    # Final write
    fem = fem_solver_combine(par, E)
    compliance = fem.compliance
    u = fem.u
    σ = fem.σ
    ε = fem.ε
    E_node = fem.E_node
    U = fem.U
    full_path = joinpath(directory, name_of_file)

    VTKGridFile(full_path, dh) do vtk
        write_solution(vtk, dh, u)
        for (j, key) in enumerate(("11", "22", "12"))
            write_cell_data(vtk, σ[j], "stress_" * key)
        end
        for (j, key) in enumerate(("11", "22", "12"))
            write_cell_data(vtk, ε[j], "strain_" * key)
        end
        write_cell_data(vtk, E, "Young's modulus")
        write_cell_data(vtk, U, "Strain Energy")
        write_node_data(vtk, E_node, "Nodal Young's modulus")
        Ferrite.write_cellset(vtk, grid)
    end

    return compliance

end
################################################################################
function top_2d(::Type{Tone} , par::DynamicParams, E, k, γ, η ,volfrac, name_of_file::String, directory::String)
    grid = par.grid
    dh = par.dh
    #E = par.E
    # nx = par.nx ; ny = par.ny ; nz = par.nz
    tnele = par.tnele
    E0 = par.E0 ; Emin = par.Emin ; Emax = par.Emax
    #k = par.k ; γ = par.γ ; volfrac = par.vf; 
    #η = par.η; 
    ρ0 = par.ρ0
    max_itr = par.max_itr ; tol = par.tol

    loop = 1

    ## Initial FEM solve
    fem = fem_solver_combine(par, E)
    compliance = fem.compliance
    H = fem.H
    W_tot = sum(fem.U)
    strain_energy_vector = [W_tot, W_tot * 10]
    A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]
    println("Iter $loop: C = $compliance, ΔR = $A")
    ## Iterative optimization loop
    while abs(A) > tol && loop <= max_itr
        fem = fem_solver_combine(par, E)
        compliance = fem.compliance
        H = fem.H
        W_tot = sum(fem.U)
        # Enew = update_young_Tone(k, E, H, Emax, Emin, E0, γ)
        # ρ = transfer_to_density(Enew, E0, ρ0, γ)
        # ρnew = filter_density_to_vf!(ρ, volfrac, tnele, η)
        # Enew_frac = transfer_to_young(ρnew, E0, ρ0, γ, Emin, Emax)

        if volfrac == 0.0
            Enew = update_young_Tone(k, E, H, Emax, Emin, E0, γ)
            Enew_frac = Enew
        elseif volfrac > 0.0
            Enew = update_young_Tone(k, E, H, Emax, Emin, E0, γ)
            ρ = transfer_to_density(Enew, E0, ρ0, γ)
            ρnew = filter_density_to_vf!(ρ, volfrac, tnele, η)
            Enew_frac = transfer_to_young(ρnew, E0, ρ0, γ, Emin, Emax)
        else
            error("Invalid value for volfrac")
        end

        # Update E in par so that fem_solver uses the updated material distribution
        #par.E = Enew_frac
        E = Enew_frac

        fem = fem_solver_combine(par, E)
        compliance = fem.compliance
        H = fem.H
        W_tot = sum(fem.U)

        strain_energy_vector[1] = strain_energy_vector[2]
        strain_energy_vector[2] = W_tot
        A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]

        loop += 1
        println("Iter $loop: C = $compliance, ΔR = $A")
    end

    ## Handle termination
    if loop > max_itr
        compliance = -1
    end

    # Final write
    fem = fem_solver_combine(par, E)
    compliance = fem.compliance
    u = fem.u
    σ = fem.σ
    ε = fem.ε
    E_node = fem.E_node
    U = fem.U
    full_path = joinpath(directory, name_of_file)

    VTKGridFile(full_path, dh) do vtk
        write_solution(vtk, dh, u)
        for (j, key) in enumerate(("11", "22", "12"))
            write_cell_data(vtk, σ[j], "stress_" * key)
        end
        for (j, key) in enumerate(("11", "22", "12"))
            write_cell_data(vtk, ε[j], "strain_" * key)
        end
        write_cell_data(vtk, E, "Young's modulus")
        write_cell_data(vtk, U, "Strain Energy")
        write_node_data(vtk, E_node, "Nodal Young's modulus")
        Ferrite.write_cellset(vtk, grid)
    end

    return compliance

end
################################################################################
function top_2d(::Type{Ttwo} , par::DynamicParams, E, B_Δt, γ, η ,volfrac, name_of_file::String, directory::String)
    grid = par.grid
    dh = par.dh
    #E = par.E
    # nx = par.nx ; ny = par.ny ; nz = par.nz
    tnele = par.tnele
    E0 = par.E0 ; Emin = par.Emin ; Emax = par.Emax
    #k = par.k ; γ = par.γ ; volfrac = par.vf; 
    #η = par.η; 
    ρ0 = par.ρ0
    max_itr = par.max_itr ; tol = par.tol

    loop = 1

    ## Initial FEM solve
    fem = fem_solver_combine(par, E)
    compliance = fem.compliance
    H = fem.H
    W_tot = sum(fem.U)
    strain_energy_vector = [W_tot, W_tot * 10]
    A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]
    println("Iter $loop: C = $compliance, ΔR = $A")
    ## Iterative optimization loop
    while abs(A) > tol && loop <= max_itr
        fem = fem_solver_combine(par, E)
        compliance = fem.compliance
        H = fem.H
        W_tot = sum(fem.U)
        # Enew = update_young_Ttwo(E, H, Emax, Emin, E0, γ, B_Δt)
        # ρ = transfer_to_density(Enew, E0, ρ0, γ)
        # ρnew = filter_density_to_vf!(ρ, volfrac, tnele, η)
        # Enew_frac = transfer_to_young(ρnew, E0, ρ0, γ, Emin, Emax)

        if volfrac == 0.0
            Enew = update_young_Ttwo(E, H, Emax, Emin, E0, γ, B_Δt)
            Enew_frac = Enew
        elseif volfrac > 0.0
            Enew = update_young_Ttwo(E, H, Emax, Emin, E0, γ, B_Δt)
            ρ = transfer_to_density(Enew, E0, ρ0, γ)
            ρnew = filter_density_to_vf!(ρ, volfrac, tnele, η)
            Enew_frac = transfer_to_young(ρnew, E0, ρ0, γ, Emin, Emax)
        else
            error("Invalid value for volfrac")
        end

        # Update E in par so that fem_solver uses the updated material distribution
        #par.E = Enew_frac
        E = Enew_frac

        fem = fem_solver_combine(par, E)
        compliance = fem.compliance
        H = fem.H
        W_tot = sum(fem.U)

        strain_energy_vector[1] = strain_energy_vector[2]
        strain_energy_vector[2] = W_tot
        A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]

        loop += 1
        println("Iter $loop: C = $compliance, ΔR = $A")

    end

    ## Handle termination
    if loop > max_itr
        compliance = -1
    end

    # Final write
    fem = fem_solver_combine(par, E)
    compliance = fem.compliance
    u = fem.u
    σ = fem.σ
    ε = fem.ε
    E_node = fem.E_node
    U = fem.U
    full_path = joinpath(directory, name_of_file)

    VTKGridFile(full_path, dh) do vtk
        write_solution(vtk, dh, u)
        for (j, key) in enumerate(("11", "22", "12"))
            write_cell_data(vtk, σ[j], "stress_" * key)
        end
        for (j, key) in enumerate(("11", "22", "12"))
            write_cell_data(vtk, ε[j], "strain_" * key)
        end
        write_cell_data(vtk, E, "Young's modulus")
        write_cell_data(vtk, U, "Strain Energy")
        write_node_data(vtk, E_node, "Nodal Young's modulus")
        Ferrite.write_cellset(vtk, grid)
    end

    return compliance

end
