### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ aa1c17e2-d62c-11ed-0431-810500d65ead
begin
	using DrWatson
	@quickactivate :SAME
	using SAME
end

# ╔═╡ a05bab50-065d-4bec-a9a7-7ca0c60db30f
begin
	using GraphPlot
	using PlutoUI
	using Plots
end

# ╔═╡ 2123deb2-8321-4e80-b1da-9f04157cb945
begin
N_individuals = 15
N_groups = 2  

p_dist = Dirac(N_individuals)                                              
g_dist = Dirac(0)

# edge_list,vertex_attribues = SAME.generate_hypergraph_bipartite_edge_list(N_groups,N_individuals,p_dist,g_dist)
# graph,vertex_attribues = SAME.make_hypergraph(N_groups,N_individuals,p_dist,g_dist) 
 graph,hypergraph_meta = SAME.make_SBM(N_groups,N_individuals,p_dist,g_dist,0.00)
end

# ╔═╡ 79e30ac9-db2a-47ef-8d66-7e45288c9e95
md"""
### Hello World
Lets take a look at the dynamics of the voter model
"""

# ╔═╡ 4384683d-63be-45c5-8e05-7815b22a95cb
begin 
sn = SAME.init_spin_network(graph)
my_vm = SAME.VoterModel(sn,SAME.voter_model_update)
vm_ts = SAME.run_model(my_vm)
end

# ╔═╡ f5d26520-e1fa-4ee0-9cbb-264256682b53


# ╔═╡ 5ddb1224-3260-4ed6-abfd-97bc0f008341


# ╔═╡ f4abc3d4-ebdd-47c1-8247-bd703b36d4df
md"""
# Voter Model Exploration

What do the dynamics look like for the standard voter model? Lets look at the graph and the average magnetization over time? 
$(@bind a Slider(1:1000))
"""

# ╔═╡ 938ce8f6-4938-4402-bdcf-b9db640adc46
begin
	my_sn = vm_ts[1].spin_network
	my_spin_vals = my_sn.spin_vals
	my_graph = my_sn.g
	nodecolordict = Dict(-1=>"lightseagreen", 1 => "orange")
	nodefillc = map(x ->nodecolordict[x],my_spin_vals)
	graph_plot = gplot(my_graph,nodefillc = nodefillc,layout = circular_layout)
end

# ╔═╡ 5db2496f-d819-4d54-94c4-16069af07e54
md """Magnetization over time"""

# ╔═╡ 969d2128-33c3-469a-8b3c-d3c79394a892
begin
	vm_mag = map(x->SAME.get_magnetization(x),vm_ts)
	mag_plot = plot(1:a,vm_mag[1:a])
end

# ╔═╡ edb6d5d3-047b-4055-8f02-6bdbc7974ecc
md"""
# Ising Model
What do the dynamics look like for the standard ising model? Lets look at the graph and the average magnetization over time? 
$(@bind ising_t Slider(1:1000))
"""


# ╔═╡ 078800ba-5fd4-4baa-bdc6-765ff72ebeb5
begin 
my_lattice = SAME.grid((100,100))
ising_sn = SAME.init_spin_network(my_lattice)
my_im = SAME.IsingModel(ising_sn,SAME.metropolis_hastings_update, 0.0,0.1,1.0)
im_ts = SAME.run_model(my_im)
end

# ╔═╡ 80e7b1d6-c39b-4faf-9f9d-1a7aeb26d099
begin
	my_im_sn = im_ts[ising_t].spin_network
	my_im_spin_vals = my_im_sn.spin_vals
	my_im_graph = my_im_sn.g
	im_nodefillc = map(x ->nodecolordict[x],my_im_spin_vals)
	gplot(my_im_graph,nodefillc = nodefillc)
end

# ╔═╡ ff68f5ae-9bbe-41e0-ad66-82f3ad39b9cb
begin
	im_mag = map(x->SAME.get_magnetization(x),im_ts)
	plot(1:ising_t,im_mag[1:ising_t])
end

# ╔═╡ 0e585145-2b66-4233-a473-872d2fa7aa5d
SAME.grid()

# ╔═╡ Cell order:
# ╠═aa1c17e2-d62c-11ed-0431-810500d65ead
# ╟─a05bab50-065d-4bec-a9a7-7ca0c60db30f
# ╠═2123deb2-8321-4e80-b1da-9f04157cb945
# ╠═79e30ac9-db2a-47ef-8d66-7e45288c9e95
# ╠═4384683d-63be-45c5-8e05-7815b22a95cb
# ╠═f5d26520-e1fa-4ee0-9cbb-264256682b53
# ╠═5ddb1224-3260-4ed6-abfd-97bc0f008341
# ╠═f4abc3d4-ebdd-47c1-8247-bd703b36d4df
# ╠═938ce8f6-4938-4402-bdcf-b9db640adc46
# ╠═5db2496f-d819-4d54-94c4-16069af07e54
# ╠═969d2128-33c3-469a-8b3c-d3c79394a892
# ╠═edb6d5d3-047b-4055-8f02-6bdbc7974ecc
# ╠═078800ba-5fd4-4baa-bdc6-765ff72ebeb5
# ╠═80e7b1d6-c39b-4faf-9f9d-1a7aeb26d099
# ╠═ff68f5ae-9bbe-41e0-ad66-82f3ad39b9cb
# ╠═0e585145-2b66-4233-a473-872d2fa7aa5d
