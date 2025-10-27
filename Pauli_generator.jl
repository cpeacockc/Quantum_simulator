#Code to take in strings and generate pauli-based Hamiltonians
using SparseArrays, LinearAlgebra

function pauli_expand(paulistr::Union{Vector{Int64},Vector{Float64}},flag::String)
    I=spdiagm([1,1])
    local X
    local Y
    local Z
    if flag=="spin"
        X=Complex.(sparse([2,1],[1,2],[1/2,1/2]))
        Y=Complex.(sparse([2,1],[1,2],[im/2,-im/2]))
        Z=Complex.(spdiagm([1/2,-1/2])) 
    elseif flag=="pauli"
        X=Complex.(sparse([2,1],[1,2],[1,1]))
        Y=Complex.(sparse([2,1],[1,2],[im,-im]))
        Z=Complex.(spdiagm([1,-1]))
    else
        println("invalid spinflag (options: spin, pauli)") 
    end
    P = X+im*Y
    M = X-im*Y
    N = length(paulistr)
    BigOp=0
    for i in 1:1:N
        n=paulistr[i]
        if n==0
            Op=I
        elseif n==1
            Op=X
        elseif n==2
            Op=Y
        elseif n==3
            Op=Z
        elseif n==4
            Op=P
        elseif n==5
            Op=M
        end

        if i==1
            BigOp=copy(Op)
        else
            BigOp = kron(Op,BigOp)
        end
    end

    #if real(BigOp)==BigOp
    #    BigOp = real(BigOp)
    #end

    return BigOp
end

function Op_create(N::Int64,Op::String,site::Int64,flag::String)
    str = zeros(N)
    if Op == "X"
        str[site]=1
    elseif Op == "Y"
        str[site]=2
    elseif Op == "Z"
        str[site]=3
    elseif Op =="+"
        str[site]=4
    elseif Op == "-"
        str[site]=5
    end
    return pauli_expand(str,flag)
end



function Product_state(N::Int64,state_flag::String)
    up=[1, 0]
    dn=[0, 1]
    if state_flag=="Neel"
        state=copy(up)
        for i in 2:N
            if iseven(i)==true
                state=kron(state,dn)
            else
                state=kron(state,up)
            end
        end
    elseif state_flag=="Up"
        state=copy(up)
        for i in 2:N
            state=kron(state,up)
        end
    elseif state_flag=="Dn"
        state=copy(dn)
        for i in 2:N
            state=kron(state,dn)
        end
    elseif state_flag=="Rand"
        state = rand(2^N)
    end
    return state./norm(state) #normalize
end

#Generate simple unitary
U_t(t::Real,H::AbstractMatrix) = exp(-im .* Matrix(H) .* t)

# Simple commutator and anti-commutator functions
comm(A::AbstractArray,B::AbstractArray) = A*B - B*A
acomm(A::AbstractArray,B::AbstractArray) = A*B + B*A



