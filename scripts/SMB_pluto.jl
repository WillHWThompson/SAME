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

# ╔═╡ f486d60e-d64a-11ed-2d05-cd53ca1d7bba
begin
using DrWatson
@quickactivate :SAME
using SAME 
using GraphPlot
using Plots
using PlutoUI
using OrdinaryDiffEq, Plots, RecursiveArrayTools
using OrdinaryDiffEq:ODESolution
end

# ╔═╡ 68633cd9-a93b-42fd-bee0-d99524204207


# ╔═╡ 6bd4b430-48fd-4ac8-9968-0e22999f9d51
begin
N_individuals = 15
N_groups = 3
occ_prob = 0.0
end

# ╔═╡ aa72a96b-24f3-4867-8e36-934fe3128733
begin
p_dist = Dirac(N_individuals)                                              
g_dist = Dirac(0)
end

# ╔═╡ 6bd241ab-0368-4949-b759-7e1aba8ea01c
begin
p_range = 1:10:500
a = p_range
p_plot = plot(p_range,pdf.(p_dist,p_range),title = "Clique Size Distribution")
g_range = 1:20
g_plot = plot(g_range,pdf.(g_dist,g_range),title = "Group Membership Distribution")

plot(p_plot,g_plot,share_x_axis = false)
end

# ╔═╡ 221ee34c-ab6d-4513-9373-9d5b75416375
md"""
### HyperGraph Structure
"""

# ╔═╡ a66983c7-b296-49ac-85c9-3dad07ca8985
begin
# edge_list,vertex_attribues = SAME.generate_hypergraph_bipartite_edge_list(N_groups,N_individuals,p_dist,g_dist)
# graph,vertex_attribues = SAME.make_hypergraph(N_groups,N_individuals,p_dist,g_dist)
unipartite,hypergraph_meta = #SAME.make_SBM(N_groups,N_individuals,p_dist,g_dist,occ_prob)

unipartite,hypergraph_meta =SAME.make_unipartite_community_graph(N_groups,N_individuals,p_dist,g_dist)

	
nodecolordict = Dict(2=>"lightseagreen", 1 => "orange")
attribs_list = sort(hypergraph_meta.vertex_attributes) |> values |> collect
nodefillc = map(x ->nodecolordict[x],attribs_list)
gplot(hypergraph_meta.hypergraph,nodefillc = nodefillc)
end

# ╔═╡ 239a9d79-ce90-4f07-b0f0-655d4a790afd
hypergraph_meta

# ╔═╡ fd49ec62-6c09-4ad5-a92c-309140693c93
md"""
### Unipartite Projection
"""

# ╔═╡ 8277233b-a9d9-40d5-9e2c-0f470c884df0
begin
	layout=(args...)->spring_layout(args...; C=20)
	gplot(unipartite,nodesize = [0.1 for i in 1:nv(unipartite)],layout = layout)
end

# ╔═╡ d3c6d3ba-ef2c-4fc9-b2ed-6f9232af75a3


# ╔═╡ ca91a4bf-5b6a-4ffa-b2b6-4c7a42a2e4aa
begin
	print("hello")
	G_0 = SAME.initialize_u0(N_individuals,p_dist)
	dG = zeros(N_individuals,N_individuals)
	SAME.voter_model_ame!(dG,G_0,0,1)
end

# ╔═╡ 62c2e252-3166-40f5-90ef-cffefe7bed1e
begin
	t_max = 100
	p = [0]
	u₀ = SAME.initialize_u0(N_individuals,p_dist)
	t_span = (0.,t_max)
	prob = ODEProblem(SAME.voter_model_ame!,u₀,t_span,p)
	sol = solve(prob,Tsit5(),saveat=1)
end

# ╔═╡ ed6985da-7c2f-461d-add2-ab5e40eb07ab
@bind time Slider(1:t_max)

# ╔═╡ 0dd823bb-d3e2-4f6f-adc0-297e00540c98
begin
	heatmap(sol[time],xlabel = "Number of Ups",ylabel = "GroupSize")
	plot!(1:N_individuals,1:N_individuals)
end

# ╔═╡ eda3e3d4-e7df-4736-9a04-3aa77973d1a8
plot(map(x ->sum(x),sol))

# ╔═╡ Cell order:
# ╠═f486d60e-d64a-11ed-2d05-cd53ca1d7bba
# ╠═68633cd9-a93b-42fd-bee0-d99524204207
# ╠═6bd4b430-48fd-4ac8-9968-0e22999f9d51
# ╠═aa72a96b-24f3-4867-8e36-934fe3128733
# ╟─6bd241ab-0368-4949-b759-7e1aba8ea01c
# ╟─221ee34c-ab6d-4513-9373-9d5b75416375
# ╠═239a9d79-ce90-4f07-b0f0-655d4a790afd
# ╠═a66983c7-b296-49ac-85c9-3dad07ca8985
# ╟─fd49ec62-6c09-4ad5-a92c-309140693c93
# ╠═8277233b-a9d9-40d5-9e2c-0f470c884df0
# ╠═d3c6d3ba-ef2c-4fc9-b2ed-6f9232af75a3
# ╠═ca91a4bf-5b6a-4ffa-b2b6-4c7a42a2e4aa
# ╠═62c2e252-3166-40f5-90ef-cffefe7bed1e
# ╠═0dd823bb-d3e2-4f6f-adc0-297e00540c98
# ╠═ed6985da-7c2f-461d-add2-ab5e40eb07ab
# ╠═eda3e3d4-e7df-4736-9a04-3aa77973d1a8
