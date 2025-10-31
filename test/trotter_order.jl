using Revise
using BenchmarkTools
using LinearAlgebra

using QuantumChannelTrajectories

Nx = 4
Ny = 4

N = Nx * Ny

for staggered in (true, false)
    for alternating in (true, false)
        order = default_circuit_order(Nx, Ny, staggered, alternating)
        valid_circuit, message = QuantumChannelTrajectories._check_trotter_order(Nx, Ny, order)
        println(message)
        @assert valid_circuit "Generated Trotter order is invalid for staggered=$staggered, alternating=$alternating"

        println("Trotter order for staggered=$staggered, alternating=$alternating is:")
        println(order)
        println("---------------------------------------------------")
    end
end