using LinearAlgebra

# Function for diagonalizing H and time evolving a state, measuring ops
function Diag_TE_measure(psi::AbstractVector, H::AbstractMatrix, t_array::Vector{<:Real}, Ops::AbstractArray)
    eigH=eigen(H)
    Measure_t_array=Matrix{ComplexF64}(undef, length(t_array), length(Ops))

    c_n = eigH.vectors' * psi  # Project psi onto eigenbasis

    for (ti,t) in (enumerate(t_array))
        phase_factors = exp.(-im * t .* eigH.values)
        psi_t = eigH.vectors * (phase_factors .* c_n)

        Measure_t_array[ti,:] = ([(psi_t'*Ops[i]*psi_t) for i in eachindex(Ops)])

    end
    return Measure_t_array
end

#Function for time evolving a wave function from diagonalizing H
function Diag_TE(psi::AbstractVector, H::AbstractMatrix, t_array::Vector{<:Real})
    eigH = eigen(H)
    c_n = eigH.vectors' * psi  # Project psi onto eigenbasis
    psi_t_array = Matrix{ComplexF64}(undef, length(t_array), length(psi))

    for (ti, t) in (enumerate(t_array))
        phase_factors = exp.(-im * t .* eigH.values)
        psi_t = eigH.vectors * (phase_factors .* c_n)
        psi_t_array[ti, :] = psi_t
    end

    return psi_t_array
end
