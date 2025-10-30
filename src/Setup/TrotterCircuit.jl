export create_circuit
export create_circuit_Adam
export circuit_order


""" Checks if the tuples in each layer are disjoint pairs of site indices. I.e., no site index appears more than once in a layer.
    Also checks that the site indices are within the valid range [1, Nx*Ny], and that the tuples are nearest-neighbours.
    
    Returns true if the order is valid, false otherwise.
"""
function check_trotter_order(Nx::Int, Ny::Int, order)::Bool
    
    for layer in order
        seen_indices = Set{Int}()
        for (i, j) in layer
            if i in seen_indices || j in seen_indices
                return false, "Overlapping indices in layer detected."
            end
            push!(seen_indices, i)
            push!(seen_indices, j)
        end
    end

    N = Nx * Ny
    for layer in order
        for (i, j) in layer
            if i < 1 || i > N || j < 1 || j > N
                return false, "Indices ($i, $j) out of bounds."  # Indices out of bounds
            end
            if abs(i - j) != 1 && abs(i - j) != Nx
                return false, "Indices ($i, $j) are not nearest-neighbours."  # Not nearest neighbours
            end
        end
    end

    return true, message

end


"""
    create_circuit(Nx::Int, Ny::Int, Order)

    Creates a Trotter circuit for a 2D system with dimensions Nx and Ny,
    following the specified order of two-site interactions.

    Parameters:
    - `Nx::Int`: Number of sites in the x-direction.
    - `Ny::Int`: Number of sites in the y-direction.
    - `Order::Vector{Vector{Tuple{Int, Int}}}`: A vector of layers, where each layer is a vector of tuples representing pairs of site indices to interact.
    - `B::Float64`: Magnetic field strength (default is 0.0).
    - `V::Float64`: Interaction strength (default is 0.0).
    - `fermions::Bool`: Whether to use fermionic statistics (default is false).
"""
function create_circuit(Nx::Int, Ny::Int, order; B::Float64 = 0.0, V::Float64 = 0.0, fermions::Bool = false)
    
    trotter_valid, message = check_trotter_order(Nx, Ny, order)
    if !trotter_valid
        error("Invalid Trotter order: $message")
    end

    N::Int = Nx * Ny

    num_layers = length(order)
    circuit_layers = [spzeros(Complex{Float64}, 2^N, 2^N) for _ in 1:num_layers]  # Create a Hamiltonian for each layer
    
    fill_operator = fermions ? PauliZ : sparse(I, 2, 2)  # Operator to fill the JW string with (fermions vs bosons)
    
    for layer in 1:num_layers
        num_indices = length(order[layer])
        for indices in 1:num_indices
            x1, x2 = order[layer][indices]
            nx::Int = mod(x1-1,Nx)+1
            ny::Int = div(x1-1,Nx)
            
            if abs(x1-x2)==1

                local_operator = -exp(im*B*ny)*kron(Sigma_plus,Sigma_minus) - exp(-im*B*ny)*kron(Sigma_minus,Sigma_plus) + V*kron(density_operator,density_operator)

                row_operator = kron( 
                    sparse(I, 2^(nx-1), 2^(nx-1)), 
                    local_operator, 
                    sparse(I, 2^(Nx-nx-1), 2^(Nx-nx-1)) 
                )

                circuit_layers[layer] += kron( 
                    sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))), 
                    row_operator, 
                    sparse(I, 2^(N - Nx - Nx*(ny-1)), 2^(N - Nx - Nx*(ny-1)))
                )        
            
            elseif  abs(x1-x2)==Nx

                local_operator = V*kron(density_operator, sparse(I,2^(Nx-1),2^(Nx-1)), density_operator)  # build JW
                local_operator -= kron(Sigma_plus, fill(fill_operator,Nx-1)... ,Sigma_minus) 
                local_operator -= kron(Sigma_minus, fill(fill_operator,Nx-1)..., Sigma_plus) 

                row_operator = kron( 
                    sparse(I, 2^(nx-1), 2^(nx-1)), 
                    local_operator, 
                    sparse(I, 2^(Nx-nx), 2^(Nx-nx)) 
                )

                circuit_layers[layer] += kron(
                    sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))), 
                    row_operator, 
                    sparse(I, 2^(N - 2*Nx - Nx*(ny-1)), 2^(N - 2*Nx - Nx*(ny-1))) 
                ) 
                          
            else
                error("Indices ($x1, $x2) are not nearest-neighbours. This wasn't caught in the initial check.")
            end
        end
    end

    # manual cleanup
    row_operator = nothing
    local_operator = nothing

    return circuit_layers
end


function default_circuit_order(Nx::Int, Ny::Int, staggered::Bool, alternating::Bool)

    layers=[[] for _ in 1:4]

    # horizontal bonds
    for ny in 1:Ny
        for nx in 1:Nx-1
            index = nx + (ny-1)*Nx
            
            ## NEEDS FINISHING ###

        end
    end

end


function circuit_order(Nx::Int,Ny::Int; type="Anna")
    
    


    
    if type=="Anna"
        for row in 1:Ny
            for col in 1:Nx
                index = col+(row-1)*Nx
                col<Nx && (col+row)%2==0 ? push!(layers[1], (index,index+1)) : nothing    
                row<Ny && (col+row)%2==0 ? push!(layers[3], (index,index+Nx)) : nothing    
                col<Nx && (col+row)%2==1 ? push!(layers[2], (index,index+1)) : nothing    
                row<Ny && (col+row)%2==1 ? push!(layers[4], (index,index+Nx)) : nothing     
            end
        end
    elseif type=="Adam"
        for row in 1:Ny
            for col in 1:Nx
                index = col+(row-1)*Nx
                col<Nx && col%2==1 ? push!(layers[1], (index,index+1)) : nothing    
                col<Nx && col%2==0 ? push!(layers[3], (index,index+1)) : nothing    
                row<Ny && row%2==1 ? push!(layers[2], (index,index+Nx)) : nothing    
                row<Ny && row%2==0 ? push!(layers[4], (index,index+Nx)) : nothing     
            end
        end
    else   
        for row in 1:Ny
            for col in 1:Nx
                index = col+(row-1)*Nx
                col<Nx && col%2==1 ? push!(layers[1], (index,index+1)) : nothing    
                col<Nx && col%2==0 ? push!(layers[2], (index,index+1)) : nothing    
                row<Ny && row%2==1 ? push!(layers[3], (index,index+Nx)) : nothing    
                row<Ny && row%2==0 ? push!(layers[4], (index,index+Nx)) : nothing     
            end
        end
    end
    
    return layers
end

