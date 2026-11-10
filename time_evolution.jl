using ITensors, ITensorMPS
using LinearAlgebra
using JLD2
using SHA

include(joinpath(@__DIR__, "DMRG.jl"))

function write_correlations(io, step, time, corr)
    for site1 in axes(corr, 1), site2 in axes(corr, 2)
        println(io, step, ",", time, ",", site1, ",", site2, ",", corr[site1, site2])
    end
    flush(io)
end

function write_response_correlations(io, step, time, excitation_site, correlations)
    for (site, correlation) in enumerate(correlations)
        println(io, step, ",", time, ",", excitation_site, ",", site, ",", correlation)
    end
    flush(io)
end

function write_magnetization(io, step, time, magz)
    for site in eachindex(magz)
        println(io, step, ",", time, ",", site, ",", magz[site])
    end
    flush(io)
end

function local_response_correlations(psi0, phi, sites, excitation_norm, energy0, time)
    amplitude = excitation_norm^2 * exp(1im * energy0 * time)
    return [
        amplitude * inner(psi0, apply(op("S+", sites, site), phi; cutoff=1e-12))
        for site in eachindex(sites)
    ]
end

function save_time_checkpoint(
    checkpoint_file,
    phi,
    elapsed,
    step,
    step_dt,
    requested_dt,
    excitation_site,
    excitation_norm,
    energy0,
    bond_digest,
)
    save(
        checkpoint_file,
        "phi", phi,
        "elapsed", elapsed,
        "step", step,
        "step_dt", step_dt,
        "requested_dt", requested_dt,
        "excitation_site", excitation_site,
        "excitation_norm", excitation_norm,
        "energy0", energy0,
        "bond_digest", bond_digest,
    )
end

function run_time_evolution(
    bond_file::String=BOND_FILE;
    total_time::Float64=1.0,
    dt::Float64=0.01,
    excitation_site::Union{Nothing,Int}=nothing,
    outdir::String=joinpath(DMRG_OUTPUT_DIR, "time_evolution"),
    resume::Bool=false,
)
    isfinite(total_time) && total_time > 0 ||
        throw(ArgumentError("total_time must be finite and positive"))
    isfinite(dt) && dt > 0 ||
        throw(ArgumentError("dt must be finite and positive"))

    bonds = read_bonds(bond_file)
    psi_file = joinpath(DMRG_OUTPUT_DIR, "psi0.jld2")
    psi0, energy0 = load(psi_file, "psi0", "energy")
    sites = siteinds(psi0)
    N = length(sites)
    all(max(bond.site1, bond.site2) <= N for bond in bonds) ||
        throw(ArgumentError("bonds in $bond_file reference sites outside the ground state"))
    H = build_hamiltonian(bonds, sites)

    j0 = isnothing(excitation_site) ? cld(N, 2) : excitation_site
    1 <= j0 <= N || throw(ArgumentError("excitation_site must be between 1 and $N"))

    phi = apply(op("S-", sites, j0), psi0; cutoff=1e-10)
    excitation_norm = norm(phi)
    isfinite(excitation_norm) && excitation_norm > 0 ||
        throw(ArgumentError("S- at site $j0 produces a zero or invalid state"))
    normalize!(phi)

    mkpath(outdir)
    corr_zz_path = joinpath(outdir, "corr_zz.csv")
    response_pm_path = joinpath(outdir, "response_pm.csv")
    mag_z_path = joinpath(outdir, "mag_z.csv")
    checkpoint_file = joinpath(outdir, "final_state.jld2")
    bond_digest = bytes2hex(sha256(read(bond_file)))

    if resume
        isfile(checkpoint_file) ||
            throw(ArgumentError("Cannot resume: checkpoint not found at $checkpoint_file"))
        checkpoint = load(checkpoint_file)
        checkpoint["bond_digest"] == bond_digest ||
            throw(ArgumentError("Cannot resume: bonds differ from the saved checkpoint"))
        isapprox(checkpoint["requested_dt"], dt; atol=1e-12, rtol=1e-12) ||
            throw(ArgumentError("Cannot resume: dt differs from the saved checkpoint"))
        isnothing(excitation_site) || excitation_site == checkpoint["excitation_site"] ||
            throw(ArgumentError("Cannot resume: excitation site differs from the checkpoint"))

        phi = checkpoint["phi"]
        elapsed = checkpoint["elapsed"]
        step = checkpoint["step"]
        step_dt = checkpoint["step_dt"]
        j0 = checkpoint["excitation_site"]
        excitation_norm = checkpoint["excitation_norm"]
        energy0 = checkpoint["energy0"]
        sites = siteinds(phi)
        N = length(sites)

        total_time > elapsed ||
            throw(ArgumentError("total_time must exceed checkpoint time $elapsed to continue"))
        additional_steps = (total_time - elapsed) / step_dt
        isapprox(additional_steps, round(additional_steps); atol=1e-9, rtol=1e-9) ||
            throw(ArgumentError("total_time must align with the checkpoint step size $step_dt"))
        for path in (corr_zz_path, response_pm_path, mag_z_path)
            isfile(path) || throw(ArgumentError("Cannot resume: output file missing: $path"))
        end
        firstline(path) = open(path, "r") do file
            readline(file)
        end
        firstline(corr_zz_path) == "step,time,site1,site2,value" ||
            throw(ArgumentError("Cannot resume: unexpected header in $corr_zz_path"))
        firstline(response_pm_path) == "step,time,source_site,site,value" ||
            throw(ArgumentError("Cannot resume: unexpected header in $response_pm_path"))
        firstline(mag_z_path) == "step,time,site,value" ||
            throw(ArgumentError("Cannot resume: unexpected header in $mag_z_path"))
        mode = "a"
    else
        nsteps = ceil(Int, total_time / dt)
        step_dt = total_time / nsteps
        elapsed = 0.0
        step = 0
        phi = apply(op("S-", sites, j0), psi0; cutoff=1e-10)
        excitation_norm = norm(phi)
        isfinite(excitation_norm) && excitation_norm > 0 ||
            throw(ArgumentError("S- at site $j0 produces a zero or invalid state"))
        normalize!(phi)
        mode = "w"
    end

    open(corr_zz_path, mode) do corr_zz_io
        open(response_pm_path, mode) do response_pm_io
            open(mag_z_path, mode) do mag_z_io
                if !resume
                    println(corr_zz_io, "step,time,site1,site2,value")
                    println(response_pm_io, "step,time,source_site,site,value")
                    println(mag_z_io, "step,time,site,value")

                    initial_corr_pm = local_response_correlations(
                        psi0, phi, sites, excitation_norm, energy0, 0.0,
                    )
                    write_response_correlations(response_pm_io, 0, 0.0, j0, initial_corr_pm)
                end

                target_step = round(Int, total_time / step_dt)
                while elapsed < total_time - 1e-12 * step_dt
                    phi = tdvp(
                        H,
                        -1im * step_dt,
                        phi;
                        nsweeps=1,
                        maxdim=200,
                        cutoff=1e-8,
                        updater_kwargs=(; tol=1e-8, krylovdim=20),
                        normalize=true,
                        reverse_step=true,
                        outputlevel=0,
                    )
                    step += 1
                    elapsed = step * step_dt

                    corr_zz = correlation_matrix(phi, "Sz", "Sz")
                    corr_pm = local_response_correlations(
                        psi0, phi, sites, excitation_norm, energy0, elapsed,
                    )
                    mag_z = expect(phi, "Sz")
                    write_correlations(corr_zz_io, step, elapsed, corr_zz)
                    write_response_correlations(response_pm_io, step, elapsed, j0, corr_pm)
                    write_magnetization(mag_z_io, step, elapsed, mag_z)

                    println("Completed TDVP step $step/$target_step at t=$elapsed")
                end
            end
        end
    end

    save_time_checkpoint(
        checkpoint_file, phi, elapsed, step, step_dt, dt, j0,
        excitation_norm, energy0, bond_digest,
    )
    println("Saved restart checkpoint to ", checkpoint_file)

    return nothing
end

function main(args::Vector{String})
    positional = filter(!=("--resume"), args)
    length(positional) <= 3 ||
        error("Usage: julia time_evolution.jl [total_time] [dt] [excitation_site] [--resume]")
    count(==("--resume"), args) <= 1 ||
        error("--resume may only be specified once")
    total_time = length(positional) >= 1 ? parse(Float64, positional[1]) : 1.0
    dt = length(positional) >= 2 ? parse(Float64, positional[2]) : 0.01
    excitation_site = length(positional) >= 3 ? parse(Int, positional[3]) : nothing

    run_time_evolution(; total_time, dt, excitation_site, resume=("--resume" in args))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
