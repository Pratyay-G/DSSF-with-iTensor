using ITensors, ITensorMPS
using LinearAlgebra
using DelimitedFiles
using JLD2

BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)
ITensors.enable_threaded_blocksparse()

const BOND_FILE = joinpath(@__DIR__, "bonds_compound.dat")
const DMRG_OUTPUT_DIR = joinpath(@__DIR__, "results_compound")

struct Bond
    site1::Int
    site2::Int
    Jzz::Float64
    Jxx::Float64
end

function read_bonds(filename::String)
    bonds = Bond[]

    for (line_number, line) in enumerate(eachline(filename))
        fields = split(strip(line))
        isempty(fields) && continue
        startswith(first(fields), "#") && continue
        lowercase.(fields) == ["site1", "site2", "jzz", "jxx"] && continue
        length(fields) == 4 ||
            throw(ArgumentError("$filename:$line_number: expected site1 site2 Jzz Jxx"))

        site1 = tryparse(Int, fields[1])
        site2 = tryparse(Int, fields[2])
        Jzz = tryparse(Float64, fields[3])
        Jxx = tryparse(Float64, fields[4])
        (site1 !== nothing && site2 !== nothing && Jzz !== nothing && Jxx !== nothing) ||
            throw(ArgumentError("$filename:$line_number: invalid bond values"))
        site1 > 0 && site2 > 0 && site1 != site2 ||
            throw(ArgumentError("$filename:$line_number: sites must be distinct positive integers"))
        isfinite(Jzz) && isfinite(Jxx) ||
            throw(ArgumentError("$filename:$line_number: couplings must be finite"))

        push!(bonds, Bond(site1, site2, Jzz, Jxx))
    end

    isempty(bonds) && throw(ArgumentError("$filename contains no bonds"))
    return bonds
end

function build_hamiltonian(bonds::Vector{Bond}, sites)
    os = OpSum()
    for bond in bonds
        i, j = bond.site1, bond.site2
        os += (bond.Jzz, "Sz", i, "Sz", j)
        os += (bond.Jxx / 2, "S+", i, "S-", j)
        os += (bond.Jxx / 2, "S-", i, "S+", j)
    end
    return MPO(os, sites)
end

function run_dmrg(bond_file::String=BOND_FILE)
    bonds = read_bonds(bond_file)
    N = maximum(max(bond.site1, bond.site2) for bond in bonds)
    sites = siteinds("S=1/2", N)
    H = build_hamiltonian(bonds, sites)

    state = [isodd(n) ? "Up" : "Dn" for n in 1:N]
    psi = randomMPS(sites, state, 20)

    nsweeps = 20
    maxdim = [8, 8, 16, 16, 32, 32, 64, 64, 100, 100,
              150, 150, 200, 200]
    cutoff = fill(1e-10, nsweeps)
    noise = [1e-8, 1e-8, 1e-8, 1e-8, 1e-9, 1e-9, 1e-10, 1e-10,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    energy, psi0 = dmrg(H, psi; nsweeps, maxdim, cutoff, noise)

    
    mkpath(DMRG_OUTPUT_DIR)

    println("DMRG energy = ", energy)
    zzcorr = correlation_matrix(psi0, "Sz", "Sz")
    magz = expect(psi0, "Sz")
    writedlm(joinpath(DMRG_OUTPUT_DIR, "Corr_ZZ.csv"), zzcorr, ',')
    writedlm(joinpath(DMRG_OUTPUT_DIR, "Mag_Z.dat"), magz, ',')
    save(joinpath(DMRG_OUTPUT_DIR, "psi0.jld2"), "psi0", psi0, "energy", energy)

    return energy, psi0
end

function main(args::Vector{String})
    isempty(args) || error("Usage: julia DMRG.jl")
    run_dmrg()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
