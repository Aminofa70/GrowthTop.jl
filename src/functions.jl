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


"""
update_young_modulus("UPM", E, H, Emax, Emin, E0, γ, k=2.0)
update_young_modulus("SIMP", E, H, Emax, Emin, E0, γ)
update_young_modulus("T1", E, H, Emax, Emin, E0, γ, k=1.5)
update_young_modulus("T2", E, H, Emax, Emin, E0, γ, B_Δt=0.1)

"""
function update_young_modulus(
    method::String, 
    E::AbstractVector{<:Real}, 
    H::AbstractVector{<:Real}, 
    Emax::Real, 
    Emin::Real, 
    E0::Real, 
    γ::Real; 
    k::Real=1.0, 
    B_Δt::Real=1.0
)
    Enew = similar(E)
    α = similar(E)
    H_mean = mean(H)
    H_std = std(H)

    for i in eachindex(E)
        part1 = H[i] * (E[i] / E0)^((γ - 1) / γ)

        if method == "UPM"
            if part1 - H_mean < 0
                α[i] = -((abs(part1 - H_mean)) / (k * H_std))^γ
            else
                α[i] = ((part1 - H_mean) / (k * H_std))^γ
            end
            Enew[i] = E[i] * (1 + α[i])
        
        elseif method == "SIMP"
            α[i] = (H[i] * (E[i] / E0)^((γ - 1) / γ))^γ
            Enew[i] = E[i] * (1 + α[i])

        elseif method == "T1"
            if part1 - H_mean < 0
                α[i] = -((abs(part1 - H_mean)) / k)^γ
            else
                α[i] = ((part1 - H_mean) / k)^γ
            end
            Enew[i] = E[i] * (1 + α[i])

        elseif method == "T2"
            H_mean = 0.0
            ρ0 = 1.0
            residual = (E0 * B_Δt * γ) / (ρ0^2)
            if part1 - H_mean < 0
                α[i] = -((abs(part1 - H_mean)))^γ
            else
                α[i] = ((part1 - H_mean))^γ
            end
            α_scalar = α[i] * (residual^γ)
            Enew[i] = E[i] + E0 * α_scalar

        else
            error("Unsupported method: $method")
        end

        Enew[i] = clamp(Enew[i], Emin, Emax)
    end

    return Enew
end

function top_2d(method::Symbol, par::DynamicParams, E, args...; volfrac, η, name_of_file, directory)
    grid, dh = par.grid, par.dh
    tnele, E0, Emin, Emax, ρ0 = par.tnele, par.E0, par.Emin, par.Emax, par.ρ0
    max_itr, tol = par.max_itr, par.tol

    loop = 1
    strain_energy_vector = [1.0, 10.0]

    while loop <= max_itr
        fem = fem_solver(par, E)
        W_tot = sum(fem.U)
        H = fem.H

        # Update modulus
        if method == :UPM
            Enew = update_young_modulus("UPM", E, H, Emax, Emin, E0, args[2]; k=args[1])
        elseif method == :SIMP
            Enew = update_young_modulus("SIMP", E, H, Emax, Emin, E0, args[1])
        elseif method == :Tone
            Enew = update_young_modulus("T1", E, H, Emax, Emin, E0, args[2]; k=args[1])
        elseif method == :Ttwo
            Enew = update_young_modulus("T2", E, H, Emax, Emin, E0, args[1]; B_Δt=args[2])
        else
            error("Unknown method $method")
        end

        Enew_frac = if volfrac == 0.0
            Enew
        else
            ρ = transfer_to_density(Enew, E0, ρ0, args[end])
            ρnew = filter_density_to_vf!(ρ, volfrac, tnele, η)
            transfer_to_young(ρnew, E0, ρ0, args[end], Emin, Emax)
        end

        strain_energy_vector = [strain_energy_vector[2], sum(fem_solver(par, Enew_frac).U)]
        A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]
        E = Enew_frac

        if abs(A) <= tol
            break
        end
        loop += 1
    end

    # Final FEM solve for output
    fem = fem_solver(par, E)
    full_path = joinpath(directory, name_of_file)
    VTKGridFile(full_path, dh) do vtk
        write_solution(vtk, dh, fem.u)
        for (j, key) in enumerate(("11", "22", "12"))
            write_cell_data(vtk, fem.σ[j], "stress_" * key)
            write_cell_data(vtk, fem.ε[j], "strain_" * key)
        end
        write_cell_data(vtk, E, "Young's modulus")
        write_cell_data(vtk, fem.U, "Strain Energy")
        write_node_data(vtk, fem.E_node, "Nodal Young's modulus")
        Ferrite.write_cellset(vtk, grid)
    end

    return fem.compliance
end

function optimiser_2D_with_combinations(method::Symbol, par::DynamicParams, E, γ_vals, vf_vals, η_vals, special_vals, name_of_file::String, directory::String)
    if !isdir(directory)
        println("Creating directory: $directory")
        mkpath(directory)
    end

    log_path = joinpath(directory, "$name_of_file.txt")
    open(log_path, "w") do file
        header = if method == :Ttwo "File γ vf η B_Δt Compliance" else "File γ vf η k Compliance" end
        write(file, header * "\n")
    end

    index = 1
    for vf in vf_vals
        for γ in γ_vals
            if vf == 0.0
                η = 0.0
                for special in special_vals
                    args = (special, γ)
                    fname = "out_$(lpad(index, 4, '0')).vtu"
                    compliance = top_2d(method, par, E, args...; volfrac=vf, η=η, name_of_file=fname, directory=directory)
                    open(log_path, "a") do f
                        write(f, "$(splitext(fname)[1]) $γ nvf $η $special $compliance\n")
                    end
                    index += 1
                end
            else
                for η in η_vals
                    for special in special_vals
                        args = (special, γ)
                        fname = "out_$(lpad(index, 4, '0')).vtu"
                        compliance = top_2d(method, par, E, args...; volfrac=vf, η=η, name_of_file=fname, directory=directory)
                        open(log_path, "a") do f
                            write(f, "$(splitext(fname)[1]) $γ $vf $η $special $compliance\n")
                        end
                        index += 1
                    end
                end
            end
        end
    end
    println("Optimization completed. Results saved to $log_path")
end
