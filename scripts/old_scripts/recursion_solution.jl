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

# ╔═╡ 00695ce9-2efa-473a-a330-e019211e4133
import Pkg; Pkg.add("LaTeXStrings")

# ╔═╡ 083db3e6-06cf-11ee-0756-ff2539054bdd
begin
using DrWatson
@quickactivate :SAME
using SAME 
using GraphPlot
using Plots
using PlutoUI
import PlutoUI: combine
using LaTeXStrings
end

# ╔═╡ b6787347-fe6d-40c5-a8d7-409a26545b5d
begin

G_i(i,n,ρ; G_0 = 1) = i == 0 ? 1 : G_0 * (n*ρ)/(i*n - i^2 +n*ρ)

function get_G_dist(n,ρ)
	i_list = 0:n
	G_dist_unnormed = map(i -> G_i(i,n,ρ),i_list )
	G_dist_sum = sum(G_dist_unnormed)
	return G_dist_unnormed ./ G_dist_sum
end
end

# ╔═╡ 47a39149-9fed-42aa-b566-773de6d08114
function make_input(values::Vector)
	return combine() do Child
		inputs = [
			md"""
			$(name.label): $(
			Child(name.label,name.slider)
			)"""
			for name in values
		]

		md"""
		### Values 
		$(inputs)
		"""
	end

end

# ╔═╡ 22b80512-7d3c-4bd3-860a-936be784e133
@bind values_1 make_input([(label = "ρ",slider = Slider(0:0.1:5,default = 1)),
							(label = "n",slider = Slider(1:100,default = 20))
])

# ╔═╡ 0b40f753-7d02-4dd3-95e4-931cfe34b5f7
values_1

# ╔═╡ a91aff74-118b-4419-a4de-b6363084fb7e
begin
	my_i_list = 0:values_1[2]
	G_dist = get_G_dist(values_1[2],values_1[1])
	plot(my_i_list,G_dist,ylim = (0,1),
	title = "")
end

# ╔═╡ a3e8749f-71b0-4ceb-ae8a-27beaccddb4e
md"""

So this does not look good. A high coupling does not lead to a flat curve, it just leads to a parabola but we know that 

$\lim_{\rho \to \infty} G_0 \frac{n\rho}{in-i^2+n\rho} = 1$

So what gives? 
"""

# ╔═╡ 03b1ac0f-9d96-470b-9dc8-8fcf0606eb83
begin
	n_value_1 = 30
	i_range_1 = 0:n_value_1
	rho_range_1 = [0.01 0.1 0.5 1 5 10]
	G_dists_list_1 = map(rho_val -> get_G_dist(n_value_1,rho_val),vec(rho_range_1))
	
	labels_1 =  "ρ =" .* string.(rho_range_1)
	my_grad = cgrad(:acton,5,categorical = true)
	
	plot(i_range_1,
		G_dists_list_1,
		label = labels_1,
		ylim = (0,0.6),
		lw = 2,
		palette = my_grad,
		title = L"$G_u$ as a function of $\rho$",
		ylabel = L"G_u",
		xlabel = L"u")
savefig("G_u_rho_dist.png")
end

# ╔═╡ 4a57ce48-d84c-4e9b-a022-c9aba2cbb042
md""" 

Okay so this did not work, but what if we don't assume that $\rho_d$ and $\rho_u$ are equal

"""

# ╔═╡ e6d2e882-9ace-4aff-92fa-6b74f69933d4


# ╔═╡ Cell order:
# ╠═083db3e6-06cf-11ee-0756-ff2539054bdd
# ╠═00695ce9-2efa-473a-a330-e019211e4133
# ╠═b6787347-fe6d-40c5-a8d7-409a26545b5d
# ╠═47a39149-9fed-42aa-b566-773de6d08114
# ╠═22b80512-7d3c-4bd3-860a-936be784e133
# ╟─0b40f753-7d02-4dd3-95e4-931cfe34b5f7
# ╠═a91aff74-118b-4419-a4de-b6363084fb7e
# ╟─a3e8749f-71b0-4ceb-ae8a-27beaccddb4e
# ╠═03b1ac0f-9d96-470b-9dc8-8fcf0606eb83
# ╠═4a57ce48-d84c-4e9b-a022-c9aba2cbb042
# ╠═e6d2e882-9ace-4aff-92fa-6b74f69933d4
