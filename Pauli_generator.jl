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

function ED_exp(H,state,N,t_array)
    if N>=14
        println("System too large")
    else
        state_array = []
        
        for t in t_array
            push!(state_array,Op_t(t,H)*state)
        end
        return state_array
    end
end

comm(A,B) = A*B - B*A
acomm(A,B) = A*B + B*A


function DM_measure(rho,flag,sites,spinflag)
    N = Int(log2(size(rho)[1]))
    meas_arr = []
    for i in sites
        str = zeros(N)
        if flag == "X"
            str[i]=1
            A = pauli_expand(str,spinflag)
            push!(meas_arr,tr(rho*A))
        elseif flag == "Y"
            str[i]=2
            A = pauli_expand(str,spinflag)
            push!(meas_arr,tr(rho*A))
        elseif flag == "Z"
            str[i]=3
            A = pauli_expand(str,spinflag)
            push!(meas_arr,tr(rho*A))
        elseif flag == "J"
            str[i]=1
            X1 = pauli_expand(str,spinflag)
            str = zeros(N)
            str[i+1]=1
            X2 = pauli_expand(str,spinflag)
            str = zeros(N)
            str[i]=2
            Y1 = pauli_expand(str,spinflag)
            str = zeros(N)
            str[i+1]=2
            Y2 = pauli_expand(str,spinflag)
            J = -2 .* (Y1*X2-X1*Y2)
            push!(meas_arr,tr(rho*J))
        end
    end
    return meas_arr
end

function state_measure(state,op,sites,spinflag)
    N = Int(log2(size(state)[1]))
    meas_arr = []
    if typeof(op)==String
        for i in sites
            str = zeros(N)
            if op == "X"
                str[i]=1
                A = pauli_expand(str,spinflag)
                push!(meas_arr,state'*A*state)
            elseif op == "Y"
                str[i]=2
                A = pauli_expand(str,spinflag)
                push!(meas_arr,state'*A*state)
            elseif op == "Z"
                str[i]=3
                A = pauli_expand(str,spinflag)
                push!(meas_arr,state'*A*state)
            elseif op == "J"
                str[i]=1
                X1 = pauli_expand(str,spinflag)
                str = zeros(N)
                str[i+1]=1
                X2 = pauli_expand(str,spinflag)
                str = zeros(N)
                str[i]=2
                Y1 = pauli_expand(str,spinflag)
                str = zeros(N)
                str[i+1]=2
                Y2 = pauli_expand(str,spinflag)
                J = -2 .*(Y1*X2-X1*Y2)
                push!(meas_arr,state'*J*state)
            else
                println("Incorrect Op string. Try inputting matrix operator")
            end
        end
        return real.(meas_arr)
    else
        return real(state'*op*state)
    end

end


function bitarr_to_int(arr)
    arr = reverse(arr)
    return sum(arr .* (2 .^ collect(length(arr)-1:-1:0)))
end


U_t(t,H) = expt(-im*H*t)


function Schrodinger_evolution(psi_0,H,ttrange)

    state_array=[]
    for t in trange
        psi = U_t(t,H)*psi_0
        push!(psi,state_array)
    end
    return state_array
end





#=
function drho(H,state,N,t_array,A_array,g_array)
    
    drho = -im .* comm(rho,H)
    for i in 1:1:length(g_array)
        A = A_array[i]
        drho += g_array[i] .* (A*rho*A' - 0.5 .* acomm(A'A,rho))
    end
    return drho
end
=#