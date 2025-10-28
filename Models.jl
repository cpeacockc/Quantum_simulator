using SparseArrays


######## Random GOE matrix ##########
function H_GOE(N::Int64)
    R_ = randn(2^N,2^N)
    R = (1/sqrt(2)) .* (R_ + transpose(R_))
    return R
end

######## Random matrix of all 2 body pauli operators ######
function Random_XYZ(N::Int64)
    Field_vec=ones(3)
    Int_vec=ones(3)
    spinflag="pauli"

    Big_Op=spzeros(2^N,2^N)
        
    for i in 1:N
        for j in 1:N
            if i<j
                str = zeros(N)
                str[i]=1;str[j]=1 #X_i X_i+1
                Big_Op += Int_vec[1]*randn() .* pauli_expand(str,spinflag)


                str = zeros(N)
                str[i]=2;str[j]=2 #Y_i Y_i+1
                Big_Op += Int_vec[2]*randn() .* pauli_expand(str,spinflag)

                str = zeros(N)
                str[i]=3;str[j]=3 #Z_i Z_i+1
                Big_Op += Int_vec[3]*randn() .* pauli_expand(str,spinflag)
            end
        end

        str = zeros(N)
        str[i]=1;#X_i
        Big_Op += Field_vec[1]*randn() .* pauli_expand(str,spinflag)

        str = zeros(N)
        str[i]=3;#Z_i
        Big_Op += Field_vec[3]*randn() .* pauli_expand(str,spinflag)

    end
        
    return real(Big_Op).*( 2^N / ((N/2)*(3N+1))) #This normalization factor makes bandwidth scaling equal to GOE

end


########## General Heisenberg/ XYZ-spin chain ###########
function H_XYZ(N::Int64,J_vec::AbstractVector,X_vec::AbstractVector,Y_vec::AbstractVector,Z_vec::AbstractVector;periodic::Bool=true,flag::String="pauli")
    #Hopping strengths are J_vec=[Jx,Jy,Jz]
    #X,Y,Z Field strengths are X_vec,Y_vec,Z_vec respectively
    
    Jx=J_vec[1]
    Jy=J_vec[2]
    Jz=J_vec[3]
    Big_Op=spzeros(2^N,2^N)
    
    for i in 1:N-1
        str = zeros(N)
        str[i]=1;str[i+1]=1 #X_i X_i+1
        Big_Op += Jx .* pauli_expand(str,flag)


        str = zeros(N)
        str[i]=2;str[i+1]=2 #Y_i Y_i+1
        Big_Op += Jy .* pauli_expand(str,flag)

        str = zeros(N)
        str[i]=3;str[i+1]=3 #Z_i Z_i+1
        Big_Op += Jz .* pauli_expand(str,flag)
    end

    if periodic==true
        str = zeros(N)
        str[N]=1;str[1]=1 #X_i X_i+1
        Big_Op += Jx .* pauli_expand(str,flag)

        str = zeros(N)
        str[N]=2;str[1]=2 #Y_i Y_i+1
        Big_Op += Jx .* pauli_expand(str,flag)

        str = zeros(N)
        str[N]=3;str[1]=3 #Z_i Z_i+1
        Big_Op += Jx .* pauli_expand(str,flag)
    end
    
    
    for i in 1:N
        str = zeros(N)
        str[i]=1;#X_i
        Big_Op += X_vec[i] .* pauli_expand(str,flag)
    end

    for i in 1:N
    
        str = zeros(N)
        str[i]=2;#Y_i
        Big_Op += Y_vec[i] .* pauli_expand(str,flag)
    end
    
    for i in 1:N
        str = zeros(N)
        str[i]=3;#Z_i
        Big_Op += Z_vec[i] .* pauli_expand(str,flag)
    end
    
    return Big_Op
end

######### Disordered Heisenberg (model for MBL) #########
H_Disordered_Heisenberg(N::Int64,J::Real,W::Real) = H_XYZ(N,J*ones(3),zeros(N),zeros(N),W .* (rand(N).*2 .-1))

######### A chaotic Ising chain #########
function H_Quantum_Ising(N::Int64,X_field::Real,Z_field::Real;periodic::Bool=true,spinflag::String="pauli")

    Big_Op = spzeros(2^N,2^N)

    for i in 1:N-1
        str=zeros(N)
        str[i]=3 #Z_i
        str[i+1]=3 #Z_i+1
        Big_Op += pauli_expand(str,spinflag)
    end

    #add periodic BCs
    if periodic==true
    str=zeros(N)
    str[1]=1
    str[N]=1
    Big_Op += pauli_expand(str,spinflag)
    end

    for i in 1:N
        str=zeros(N)
        str[i] = 1 #X_i
        Big_Op += X_field * pauli_expand(str,spinflag)

        str=zeros(N)
        str[i] = 3 #Z_i
        Big_Op += Z_field * pauli_expand(str,spinflag)
    end

    return Big_Op
end

######## Various Anderson Models ###########
function H_Anderson1D(N::Integer,W::Real,t::Real;h::Union{Nothing, AbstractVector}=nothing,periodic::Bool=true)
    
    #If an array of disorder is not given, one will be provided
    if isnothing(h)
        h = (W/2) .* (rand(N).*2 .-1) #random potentials
    end

    #Default is sparse format
    H_Anderson = spdiagm(-1 => (-t)*ones(N-1), 1 => (-t)*ones(N-1), 0=>h)

    #If periodic boundary conditions (default), connect the boundaries
    if periodic
    H_Anderson[1,N]=-t
    H_Anderson[N,1]=-t
    end
    return H_Anderson
end

function H_Anderson2D(N::Integer,W::Real,t::Real;h::Union{Nothing, AbstractVector}=nothing,periodic::Bool=true)

    H = spzeros(N^2,N^2)
    Id = spdiagm(ones(N))
        
        Big_Op = H_Anderson1D(N,0,t;periodic)
        Big_Op = kron(Big_Op,Id)
        H += Big_Op

        Big_Op = Id
        Big_Op = kron(Big_Op,H_Anderson1D(N,0,t;periodic))
        H += Big_Op
        
    if isnothing(h)
        h = (W/2) .* (rand(N^2).*2 .-1) #random potentials
    end
    H += spdiagm(h)
    return H
end

function H_Anderson3D(N::Integer,W::Real,t::Real;h::Union{Nothing, AbstractVector}=nothing,periodic::Bool=true)

    H = spzeros(N^3,N^3)
    Id = spdiagm(ones(N))
        Big_Op = H_Anderson1D(N,0,t;periodic)
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,Id)
        H += Big_Op

        Big_Op = Id
        Big_Op = kron(Big_Op,H_Anderson1D(N,0,t;periodic))
        Big_Op = kron(Big_Op,Id)
        H += Big_Op

        Big_Op = Id
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,H_Anderson1D(N,0,t;periodic))
        H += Big_Op

    if isnothing(h)
        h = (W/2) .* (rand(N^3).*2 .-1) #random potentials
    end
    H += spdiagm(h)
    return H
end

function H_Anderson4D(N::Integer,W::Real,t::Real;h::Union{Nothing, AbstractVector}=nothing,periodic::Bool=true)

    H = spzeros(N^4,N^4)
    Id = spdiagm(ones(N))
    
        Big_Op = H_Anderson1D(N,0,t;periodic)
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,Id)
        H += Big_Op

        Big_Op = Id
        Big_Op = kron(Big_Op,H_Anderson1D(N,0,t;periodic))
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,Id)
        H += Big_Op

        Big_Op = Id
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,H_Anderson1D(N,0,t;periodic))
        Big_Op = kron(Big_Op,Id)
        H += Big_Op

        Big_Op = Id
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,Id)
        Big_Op = kron(Big_Op,H_Anderson1D(N,0,t;periodic))
        H += Big_Op

    if isnothing(h)
        h = (W/2) .* (rand(N^3).*2 .-1) #random potentials
    end
    H += spdiagm(h)
    return H
end

function H_Anderson(d::Integer,N::Integer,W::Real,t::Real;h::Union{Nothing, AbstractVector}=nothing,periodic::Bool=true)

    if d==1
       return H_Anderson1D(N,W,t;h,periodic)
    elseif d==2
        return H_Anderson2D(N,W,t;h,periodic)
    elseif d==3
        return H_Anderson3D(N,W,t;h,periodic)
    elseif d==4
        return H_Anderson4D(N,W,t;h,periodic)
    else
        error("Incorrect dimension (only d=1-4 are currently included)")
    end
end