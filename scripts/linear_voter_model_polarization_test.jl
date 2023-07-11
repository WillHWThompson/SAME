using DrWatson
@quickactivate :SAME
using SAME
using DataFrames
using PlotlyJS
#INIT MODEL
N_inds = 5000
N_p = 0.5
N_groups = trunc(Int,N_inds/10 )
M_p = [0.1,0.2,0.4,0.3,0.4,0.5,0.6,0.7]
N_steps = trunc(Int,N_inds*10)

params = @strdict N_inds N_p N_groups M_p N_steps
dicts = dict_list(params)
for (i, d) in enumerate(dicts)
    @show i
    f = SAME.run_model(d)
    savename_i = savename(dicts[i],"jld2")
    @show savename_i
    merge_dict = merge(d,f)#merge metadata with simulation output
    wsave(datadir("simulations", savename_i), merge_dict)
end


