"""
    write_lattice(output_file, bonds)

Write a lattice as rows of `(site1, site2, Jzz, Jxx)` bonds. For example,
`[(1, 2, 1.0, 0.5), (2, 3, 1.0, 0.5), (3, 1, 1.0, 0.5)]`
describes a three-site loop.
"""
function write_lattice(output_file::String, bonds)
    open(output_file, "w") do file
        println(file, "site1 site2 Jzz Jxx")
        for (site1, site2, Jzz, Jxx) in bonds
            println(file, site1, " ", site2, " ", Jzz, " ", Jxx)
        end
    end
end

"""
    write_positions(output_file, positions)

Write one `x y` coordinate row per site, ordered by site number. This
headerless two-column format is compatible with the Fourier-analysis code.
"""
function write_positions(output_file::String, positions)
    open(output_file, "w") do file
        for (x, y) in positions
            isfinite(x) && isfinite(y) ||
                throw(ArgumentError("Site coordinates must be finite"))
            println(file, x, " ", y)
        end
    end
end

function spin_chain_bonds(nsites::Int, Jzz::Float64, Jxx::Float64)
    nsites >= 1 || throw(ArgumentError("The number of sites must be positive."))

    return [(site1, site1 + 1, Jzz, Jxx) for site1 in 1:(nsites - 1)]
end

function spin_chain_positions(nsites::Int)
    nsites >= 1 || throw(ArgumentError("The number of sites must be positive."))

    return [(Float64(site), 0.0) for site in 1:nsites]
end

function main(args::Vector{String})
    1 <= length(args) <= 3 ||
        error("Usage: julia lattice.jl <number_of_sites> [Jzz] [Jxx]")

    nsites = parse(Int, args[1])
    Jzz = length(args) >= 2 ? parse(Float64, args[2]) : 1.0
    Jxx = length(args) >= 3 ? parse(Float64, args[3]) : 1.0
    bonds_file = joinpath(@__DIR__, "bonds.txt")
    positions_file = joinpath(@__DIR__, "sites.dat")

    write_lattice(bonds_file, spin_chain_bonds(nsites, Jzz, Jxx))
    write_positions(positions_file, spin_chain_positions(nsites))
    println("Wrote bond information to ", bonds_file)
    println("Wrote site positions to ", positions_file)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
