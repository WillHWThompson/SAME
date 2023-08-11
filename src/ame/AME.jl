my_binomial(n,u,p_dist;ϵ = 0.5) = pdf(p_dist,n)*binomial(BigInt(n),BigInt(u))*ϵ^u * (1-ϵ)^(n-u)

function initialize_u0(N_individuals::Int,p_dist)
    """
    initialize_u0 - create an initial distribution for the AME
    inputs: 
        N_individuals: the total number of individuals in your population
        p_dist: the probability distribution used to specify the clique size, from Distributions.jl
    """
    G = zeros(N_individuals,N_individuals+1)#second indicie is the number of upspins. There can be 0 upspins in a group but never zero individuals
    #@show G
    G_Inds = CartesianIndices(G)
    G = map(inds -> G[inds[1],inds[2]] = my_binomial(inds[1],inds[2],p_dist),G_Inds)#generate the correct distribution for each occupation number
    G_norm = G ./sum(G)

    G_norm = parse.(Float64,string.(G_norm))#convert from arbitrary to fixed precision floats

    return G_norm
end

function voter_model_ame_by_index(n_prime,u_prime,G,N,U)
    u = u_prime -1
    #n = n_prime - 1
    n = n_prime 
    dG_n_u = 0
    dG_n_u -= G[n_prime,u_prime]*(u-n)*(u/n)#upflip outflux
    dG_n_u -= G[n_prime,u_prime]*(u)*((n-u)/n)#downflip outflux
    u_prime < U && (dG_n_u += G[n_prime,u_prime+1]*((u+1)*((n-u-1)/n)))#inwards downflip flux
    u > 1 && (dG_n_u+= G[n_prime,u_prime-1]*((n-u+1)*((u-1)/n)))#inwards upflip flux
    return dG_n_u
end

function voter_model_ame_2(dG,G,t,p)
    N,U  = size(G)
    n_range = 2:N
    u_range = 1:U
    MyCartesianIndex = Iterators.product(n_range,u_range)
    map(index -> (dG[index[1],index[2]] = voter_model_ame_by_index(index[1],index[2],G,N,U)),MyCartesianIndex)
    return dG
end




# function ρ(G,k_ex,type)
#     """
#     ρ - code for the moment closure
#     """
#     N,U  = size(G)
#     #track the sum of the numerator and the denomiator sepratley
#     numerator_sum = 0 
#     denominator_sum = 0 

#     for n_prime=1:N,u_prime=1:U
#         u  = u_prime-1
#         #n = n_prime-1
#         n = n_prime

#         if type == :up
#             numerator_sum += (n-u)*G[n_prime,u_prime]*(u/n) 
#             denominator_sum += G[n_prime,u_prime]*(n-u)
#         elseif type == :down
#             numerator_sum += u*G[n_prime,u_prime]*((n-u)/n) 
#             denominator_sum += G[n_prime,u_prime]*(u)
#         else
#             println("type not a valid - must be either ':up' or ':down'")
#         end
#     end
#     return k_ex * (numerator_sum/denominator_sum)
# end

unzip(d::Dict) = (;(p.first => unzip(p.second) for p in d)...)#convert Dictionary to NamedTuple

function parse_vm_args(p)
    return_dict = Dict()
    map(i -> return_dict[Symbol(i[1])] = i[2],zip(keys(p),p))#create dictionary from named tuple

    N_vals = return_dict[:N]
    if length(N_vals) > 1
        N_range = range(N_vals...)
    else
        N_range = range(1:N_vals)
    end
    U_range = range(1,maximum(N_range)+1)
    return_dict[:N_range] = N_range
    return_dict[:U_range] = U_range

    #return unzip(return_dict)
    return NamedTuple{Tuple(Symbol.(keys(return_dict)))}(values(return_dict)) 
end



function ρ(G,p,type)
    """
    ρ - code for the moment closure
    """
    # if length(size(G)) == 1
    #     U = size(G)
    #     k_ex,N_max = p
    #     N_range = range(N_max,N_max)
    # elseif length(size(G)) == 2
    #     U,N = size(G)
    #     k_ex = p[0]
    #     N_range = range(1,N_max)
    # else
    #     println("State matrix G must be either 1 or 2 dimensional")
    # end

    #track the sum of the numerator and the denomiator sepratley
    numerator_sum = 0 
    denominator_sum = 0 
    for n_prime=p.N_range,u_prime=p.U_range
        u  = u_prime-1
        #n = n_prime-1
        n = n_prime

        if type == :up
            numerator_sum += (n-u)*G[n_prime,u_prime]*(u/n)
            denominator_sum += G[n_prime,u_prime]*(n-u)
        elseif type == :down
            numerator_sum += u*G[n_prime,u_prime]*((n-u)/n)
            denominator_sum += G[n_prime,u_prime]*(u)
        else
            println("type not a valid - must be either ':up' or ':down'")
        end
    end
    return p.k_ex * (numerator_sum/denominator_sum)
end



function ρ_1d(G,p,type)
    """
    ρ - code for the moment closure
    """
    # if length(size(G)) == 1
    #     U = size(G)
    #     k_ex,N_max = p
    #     N_range = range(N_max,N_max)
    # elseif length(size(G)) == 2
    #     U,N = size(G)
    #     k_ex = p[0]
    #     N_range = range(1,N_max)
    # else
    #     println("State matrix G must be either 1 or 2 dimensional")
    # end

    U = size(G)[1]
    n_prime,k_ex = p
    #track the sum of the numerator and the denomiator sepratley
    numerator_sum = 0 
    denominator_sum = 0 
    for u_prime in 1:U
        u  = u_prime-1
        #n = n_prime-1
        n = n_prime

        if type == :up
            numerator_sum += (n-u)*G[u_prime]*(u/n)
            denominator_sum += G[u_prime]*(n-u)
        elseif type == :down
            numerator_sum += u*G[u_prime]*((n-u)/n)
            denominator_sum += G[u_prime]*(u)
        else
            println("type not a valid - must be either ':up' or ':down'")
        end
    end
    return k_ex * (numerator_sum/denominator_sum)
end



function ρ_q(G,k_ex,q,type)
    """
    ρ - code for the moment closure
    """
    N,U  = size(G)
    #track the sum of the numerator and the denomiator sepratley
    numerator_sum = 0 
    denominator_sum = 0 

    for n_prime=1:N,u_prime=1:U
        u  = u_prime-1
        #n = n_prime-1
        n = n_prime

        if type == :up
            numerator_sum += (n-u)*G[n_prime,u_prime]*(u/n)^q
            denominator_sum += G[n_prime,u_prime]*(n-u)
        elseif type == :down
            numerator_sum += u*G[n_prime,u_prime]*((n-u)/n)^q
            denominator_sum += G[n_prime,u_prime]*(u)
        else
            println("type not a valid - must be either ':up' or ':down'")
        end
    end
    return k_ex * (numerator_sum/denominator_sum)
end




function voter_model_ame!(dG,G,p,t)
    ρ_u = ρ(G,p,:up)#up flip moment closure
    ρ_d = ρ(G,p,:down)#downflip moment closure

    for n_prime=p.N_range,u_prime=p.U_range
        u  = u_prime-1
        n = n_prime

        dG[n_prime,u_prime] = 0

        if u > n 
            continue
        end


        # (u < n) && (dG[n_prime,u_prime] -= G[n_prime,u_prime]*((n-u)*(u/n)))#upflip outflux
        # (u > 0) && (dG[n_prime,u_prime] -= G[n_prime,u_prime]*((u)*((n-u)/n)))#downflip outflux
        # (u < n) && (dG[n_prime,u_prime]+= G[n_prime,u_prime+1]*((u+1)*((n-u-1)/n)))#inwards downflip flux
        # (u > 0) && (dG[n_prime,u_prime]+= G[n_prime,u_prime-1]*((n-u+1)*((u-1)/n)))#inwards upflip flux


        (u < n) && (dG[n_prime,u_prime] -= G[n_prime,u_prime]*( (n-u)*((u/n)+ρ_u) ))#upflip outflux

        (u > 0) && (dG[n_prime,u_prime] -= G[n_prime,u_prime]*( u*(((n-u)/n)+ρ_d) ))#downflip outflux
        (u < n) && (dG[n_prime,u_prime]+= G[n_prime,u_prime+1]*( (u+1)*(((n-u-1)/n)+ρ_d)))#inwards downflip flux
        (u > 0) && (dG[n_prime,u_prime]+= G[n_prime,u_prime-1]*( (n-u+1)*(((u-1)/n)+ρ_u)))#inwards upflip flux
    end
  #  return dG 
end




function voter_model_ame_1d!(dG,G,p,t)
    N,k_ex = p 
    N = Int64(N)
    n_prime = N

    
    ρ_u = ρ_1d(G,p,:up)#up flip moment closure
    ρ_d = ρ_1d(G,p,:down)#downflip moment closure
    for u_prime in 1:size(G)[1]
        u  = u_prime-1
        n = n_prime

        dG[u_prime] = 0

        (u < n) && (dG[u_prime] -= G[u_prime]*( (n-u)*((u/n)+ρ_u) ))#upflip outflux
        (u > 0) && (dG[u_prime] -= G[u_prime]*( u*(((n-u)/n)+ρ_d) ))#downflip outflux
        (u < n) && (dG[u_prime]+= G[u_prime+1]*( (u+1)*(((n-u-1)/n)+ρ_d)))#inwards downflip flux
        (u > 0) && (dG[u_prime]+= G[u_prime-1]*( (n-u+1)*(((u-1)/n)+ρ_u)))#inwards upflip flux
    end
  #  return dG 
end



function q_voter_model_ame!(dG,G,p,t)
    N,U  = size(G)

    #_,k_ex = p...#unpack avg external degree as a paramter of the model
    q,k_ex = p 
 
    ρ_u = ρ_q(G,k_ex,q,:up)#up flip moment closure
    ρ_d = ρ_q(G,k_ex,q,:down)#downflip moment closure

    for n_prime=1:N,u_prime=1:U
        u  = u_prime-1
        #n = n_prime-1
        n = n_prime

        dG[n_prime,u_prime] = 0

        # if n == 0
        #     continue
        # end

        if u > n 
            continue
        end


        # (u < n) && (dG[n_prime,u_prime] -= G[n_prime,u_prime]*((n-u)*(u/n)))#upflip outflux
        # (u > 0) && (dG[n_prime,u_prime] -= G[n_prime,u_prime]*((u)*((n-u)/n)))#downflip outflux
        # (u < n) && (dG[n_prime,u_prime]+= G[n_prime,u_prime+1]*((u+1)*((n-u-1)/n)))#inwards downflip flux
        # (u > 0) && (dG[n_prime,u_prime]+= G[n_prime,u_prime-1]*((n-u+1)*((u-1)/n)))#inwards upflip flux


        (u < n) && (dG[n_prime,u_prime] -= G[n_prime,u_prime]*((n-u)* ((u/n)^q+ρ_u)) )#upflip outflux
        (u > 0) && (dG[n_prime,u_prime] -= G[n_prime,u_prime]*(((u)*(((n-u)/n)^q +ρ_d))))#downflip outflux
        (u < n) && (dG[n_prime,u_prime]+= G[n_prime,u_prime+1]*(((u+1)*(((n-u-1)/n)^q+ρ_d))))#inwards downflip flux
        (u > 0) && (dG[n_prime,u_prime]+= G[n_prime,u_prime-1]*(((n-u+1)*(((u-1)/n)^q+ρ_u))))#inwards upflip flux
    end
  #  return dG 
end