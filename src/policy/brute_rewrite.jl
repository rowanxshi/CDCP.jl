### A Pluto.jl notebook ###
# v0.16.3

using Markdown
using InteractiveUtils

# ╔═╡ f83b8457-7e46-4c6a-963f-bbc29194e528
using PlutoUI

# ╔═╡ a55af1d3-1fc5-4e49-a4f7-46e6c94a22e5
π(J, z) = begin
	C = length(J)
	var = sum(j -> J[j]*j, eachindex(J))
	z*var^1.1 - var
end

# ╔═╡ 92787176-d772-4d8d-bd59-f9dcaafb780a
function equalise_π(pair)
	J1 = pair[1]
	J2 = pair[2]
	var1 = -π(J1, 0)
	var2 = -π(J2, 0)
	(var1 - var2)/(var1^1.1 - var2^1.1)
end

# ╔═╡ 8cc94cb3-c9fd-4822-96a0-c6e51ee1bb2a
equalise_π((BitVector((true, false, true, false, true)), BitVector((false, true, false, true, false))))

# ╔═╡ 75f8c822-2d18-11ec-3de3-a3c4d35fb2bf
function simple_filter!(subinterval::Pair{NTuple{2, Float64}, Vector{BitVector}})
	(left, right) = subinterval.first
	options = subinterval.second
	n_options = length(options)
	
	max_left = maximum(J -> π(J, left), options)
	filter!(options) do option
		π(option, right) ≥ max_left
	end
	
	length(options) < n_options
end

# ╔═╡ 4938c1c0-e01c-4177-b6d4-bf385c6e5f48
first_pair(v) = (v[1], v[2])

# ╔═╡ 408eef73-a282-414d-9998-b834625661ef
function brute_branch!(subinterval::Pair{NTuple{2, Float64}, Vector{BitVector}}, z_equal)
	((left, right), options) = subinterval
	(J1, J2) = @views options[1:2]
	endpoints_left = (left, z_equal)
	endpoints_right = (z_equal, right)
	J1_right = π(J1, (right+z_equal)/2) > π(J2, (right+z_equal)/2)
	subinterval_right = endpoints_right => options
	copy_options = copy(options)
	subinterval_left = endpoints_left => copy_options

	deleteat!(options, J1_right ? 2 : 1)
	deleteat!(copy_options, J1_right ? 1 : 2)
	
	(subinterval_left, subinterval_right)
end

# ╔═╡ 00b5c17c-02ea-4c3e-b857-83f248c39802
function converge_brute_branches!((working, done)::NTuple{2, Vector{Pair{NTuple{2, Float64}, Vector{BitVector}}}})
	empty!(done)
	
	while !isempty(working)
		subinterval = last(working)
		if length(subinterval.second) == 1
			push!(done, pop!(working))
			continue
		end
	
		if !simple_filter!(subinterval)
			((left, right), options) = subinterval
			(J1, J2) = @views options[1:2]
			z_equal = @views equalise_π((J1, J2))
			if z_equal ≤ left
				deleteat!(options, (π(J1, left) > π(J2, left)) ? 2 : 1)
			elseif z_equal ≥ right
				deleteat!(options, (π(J1, right) > π(J2, right)) ? 2 : 1)
			else
				append!(working, brute_branch!(pop!(working), z_equal))
			end
		end
	end
	(; working, done)
end

# ╔═╡ 9d9af40e-1ec1-4cf9-9d37-76d6b293c870
out = @timed converge_brute_branches!(([(.50, 1.5) => [trues(5), falses(5), BitVector((true, false, true, false, true)), BitVector((false, true, false, true, false))]], Vector{Pair{NTuple{2, Float64}, Vector{BitVector}}}()))

# ╔═╡ 771b0b35-ec20-40de-82dd-0c31eea029c5
function test()
	V = falses(5)
	
	V_of_V = [V, V]
	n = @timed copy_V_of_V = copy(V_of_V)
	
	V[1] = true
	
	V_of_V[1][2] = true
	
	copy_V_of_V
	deleteat!(copy_V_of_V, 1)
	V_of_V
	n
end

# ╔═╡ c59b1031-d7b3-4021-b512-9bbe4b72c9af
test()

# ╔═╡ cc549063-50c7-41f4-878a-3cf59e85add8
@allocated (1, 2)

# ╔═╡ 64096de8-82d2-4acb-95d1-fb2ed2330020
append!([1], (1, 2))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.15"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "f6532909bf3d40b308a0f360b6a0e626c0e263a8"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.1"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "a8709b968a1ea6abc2dc1967cb1db6ac9a00dfb6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.5"

[[PlutoUI]]
deps = ["Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "633f8a37c47982bff23461db0076a33787b17ecd"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.15"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
"""

# ╔═╡ Cell order:
# ╠═f83b8457-7e46-4c6a-963f-bbc29194e528
# ╟─a55af1d3-1fc5-4e49-a4f7-46e6c94a22e5
# ╠═8cc94cb3-c9fd-4822-96a0-c6e51ee1bb2a
# ╟─92787176-d772-4d8d-bd59-f9dcaafb780a
# ╠═75f8c822-2d18-11ec-3de3-a3c4d35fb2bf
# ╠═4938c1c0-e01c-4177-b6d4-bf385c6e5f48
# ╠═408eef73-a282-414d-9998-b834625661ef
# ╠═00b5c17c-02ea-4c3e-b857-83f248c39802
# ╠═9d9af40e-1ec1-4cf9-9d37-76d6b293c870
# ╠═771b0b35-ec20-40de-82dd-0c31eea029c5
# ╠═c59b1031-d7b3-4021-b512-9bbe4b72c9af
# ╠═cc549063-50c7-41f4-878a-3cf59e85add8
# ╠═64096de8-82d2-4acb-95d1-fb2ed2330020
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
