my_binomial(n,u,p_dist;ϵ = 0.5) = pdf(p_dist,n)*binomial(Int(n),Int(u))*ϵ^u * (1-ϵ)^(n-u)

function initialize_u0(N_individuals::Int,p_dist)
    """
    initialize_u0 - create an initial distribution for the AME
    inputs: 
        N_individuals: the total number of individuals in your population
        p_dist: the probability distribution used to specify the clique size, from Distributions.jl
    """
    G = zeros(N_individuals,N_individuals)#second indicie is the number of upspins. There can be 0 upspins in a group but never zero individuals
    #@show G
    G_Inds = CartesianIndices(G)
    G = map(inds -> G[inds[1],inds[2]] = my_binomial(inds[1],inds[2],p_dist),G_Inds)#generate the correct distribution for each occupation number
    G_norm = G ./sum(G)
    return G_norm
end

function voter_model_ame_by_index(n_prime,u_prime,G,N,U)
    u = u_prime -1
    n = n_prime - 1
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




function ρ(G,k_ex,type)
    """
    ρ - code for the moment closure
    """
    N,U  = size(G)
    #track the sum of the numerator and the denomiator sepratley
    numerator_sum = 0 
    denominator_sum = 0 

    for n_prime=1:N,u_prime=1:U
        u  = u_prime-1
        n = n_prime-1

        #if n is zero - do not evaluate
        if n < 1
           continue 
        end

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
    return k_ex * (numerator_sum/denominator_sum)
end
    
function voter_model_ame!(dG,G,t,p)
    N,U  = size(G)

    k_ex = p#unpack avg external degree as a paramter of the model
 
    ρ_u = ρ(G,k_ex,:up)#up flip moment closure
    ρ_d = ρ(G,k_ex,:down)#downflip moment closure

    # @show ρ_u
    # @show ρ_d

    for n_prime=1:N,u_prime=1:U
        u  = u_prime-1
        n = n_prime-1
        dG[n_prime,u_prime] = 0

        if n == 0
            continue
        end

        if u > n 
            continue
        end


        # (u < n) && (dG[n_prime,u_prime] -= G[n_prime,u_prime]*((n-u)*(u/n)))#upflip outflux
        # (u > 0) && (dG[n_prime,u_prime] -= G[n_prime,u_prime]*((u)*((n-u)/n)))#downflip outflux
        # (u < n) && (dG[n_prime,u_prime]+= G[n_prime,u_prime+1]*((u+1)*((n-u-1)/n)))#inwards downflip flux
        # (u > 0) && (dG[n_prime,u_prime]+= G[n_prime,u_prime-1]*((n-u+1)*((u-1)/n)))#inwards upflip flux


        (u < n) && (dG[n_prime,u_prime] -= G[n_prime,u_prime]*((n-u)*(u/n)+ρ_u))#upflip outflux
        (u > 0) && (dG[n_prime,u_prime] -= G[n_prime,u_prime]*((u)*((n-u)/n)+ρ_d))#downflip outflux
        (u < n) && (dG[n_prime,u_prime]+= G[n_prime,u_prime+1]*((u+1)*((n-u-1)/n)+ρ_d))#inwards downflip flux
        (u > 0) && (dG[n_prime,u_prime]+= G[n_prime,u_prime-1]*((n-u+1)*((u-1)/n)+ρ_u))#inwards upflip flux
    end
  #  return dG 
end


