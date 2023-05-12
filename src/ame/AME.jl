my_binomial(n,u,p_dist;ϵ = 0.5) = pdf(p_dist,n)*binomial(BigInt(n),BigInt(u))*ϵ^u * (1-ϵ)^(n-u)

function initialize_u0(N_individuals::Int,p_dist)
    """
    initialize_u0 - create an initial distribution for the AME
    inputs: 
        N_individuals: the total number of individuals in your population
        p_dist: the probability distribution used to specify the clique size, from Distributions.jl
    """
    G = zeros(N_individuals,N_individuals)
    @show G
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
    
function voter_model_ame!(dG,G,t,p)
    N,U  = size(G)
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

        (u < n) && (dG[n_prime,u_prime] -= G[n_prime,u_prime]*(n-u)*(u/n))#upflip outflux
        (u > 0) && (dG[n_prime,u_prime] -= G[n_prime,u_prime]*(u)*((n-u)/n))#downflip outflux
        (u < n) && (dG[n_prime,u_prime]+= G[n_prime,u_prime+1]*((u+1)*((n-u-1)/n)))#inwards downflip flux
        (u > 0) && (dG[n_prime,u_prime]+= G[n_prime,u_prime-1]*((n-u+1)*((u-1)/n)))#inwards upflip flux
    end
    return dG 
end


