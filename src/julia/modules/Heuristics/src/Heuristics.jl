module Heuristics

using Data
using Parameters
using Statistics
using Solution
using DataStructures

export maxDegree, minDegree, dinMaxDegree, dsatur, getcoloring, remove!

function coloringHeuristics(inst::InstanceData,params::ParameterData, isConnect=false)

	(bestcor,bestcolors, bestorder) = (0,[-1 for i = 1:inst.NV], [-1 for i = 1:inst.NV])

	if params.problem != "cgcol"

		order = Heuristics.maxDegree(inst)
		println("MAX DEGREEE ORDER")
		#println(order)
		(cor, colors) = Heuristics.getcoloring(inst.inputG, order, inst.NV)
		println("NUMBER OF COLORS USED = $cor")
		ncolors_maxDegree = cor
		if cor > bestcor
			(bestcor,bestcolors, bestorder) = (cor,copy(colors), copy(order))
		end

		order = Heuristics.minDegree(inst)
		println("MIN DEGREEE ORDER")
		(cor, colors) = Heuristics.getcoloring(inst.inputG, order, inst.NV)
		println("NUMBER OF COLORS USED = $cor")
		ncolors_minDegree = cor
		if cor > bestcor
			(bestcor,bestcolors, bestorder) = (cor,copy(colors), copy(order))
		end

		order = reverse(Heuristics.minDegree(inst))
		println("SMALLEST DEGREEE LAST")
		(cor, colors) = Heuristics.getcoloring(inst.inputG, order, inst.NV)
		println("NUMBER OF COLORS USED = $cor")
		ncolors_smallestDegreeLast = cor
		if cor > bestcor
			(bestcor,bestcolors, bestorder) = (cor,copy(colors), copy(order))
		end

		order = Heuristics.dinMaxDegree(inst)
		println("DYNAMIC MAX DEGREEE ORDER")
		(cor, colors) = Heuristics.getcoloring(inst.inputG, order, inst.NV)
		println("NUMBER OF COLORS USED = $cor")
		ncolors_dinMaxDegree = cor
		if cor > bestcor
			(bestcor,bestcolors, bestorder) = (cor,copy(colors), copy(order))
		end
	end

	order = Heuristics.dsatur(inst)
	println("DSATUR ORDER")
	(cor, colors) = Heuristics.getcoloring(inst.inputG, order, inst.NV)
	println("NUMBER OF COLORS USED = $cor")
	ncolors_dSatur = cor
	if cor > bestcor 
		(bestcor,bestcolors,bestorder) = (cor,copy(colors), copy(order))
	end

	order = Heuristics.connectedMaxDegree(inst)
		println("CONNECTED MAX DEGREEE ORDER")
		(cor, colors) = Heuristics.getcoloring(inst.inputG, order, inst.NV)
		println("NUMBER OF COLORS USED = $cor")
		ncolors_connectedMaxDegree = cor
		if cor > bestcor || isConnect
			(bestcor,bestcolors,bestorder) = (cor,copy(colors),copy(order))
		end

	order = Heuristics.connectedMinDegree(inst)
		println("CONNECTED MIN DEGREEE ORDER")
		(cor, colors) = Heuristics.getcoloring(inst.inputG, order, inst.NV)
		println("NUMBER OF COLORS USED = $cor")
		ncolors_connectedMinDegree = cor
		if cor > bestcor
			(bestcor,bestcolors,bestorder) = (cor,copy(colors),copy(order))
		end

	println("Best (largest) number of colors = ",bestcor)

	open("resultsHeuristics_c.txt","a") do f
		write(f,"$(params.instName)")
	end

	if params.problem == "cgcol"
		open("resultsHeuristics_c.txt","a") do f
			write(f,";$(params.problem);$(params.method);$(params.form);$ncolors_dSatur;$ncolors_connectedMinDegree,$ncolors_connectedMaxDegree\n")
		end
	else
		open("results_h.txt","a") do f
			write(f,";$(params.problem);$(params.method);$(params.form);$ncolors_maxDegree;$ncolors_minDegree;$ncolors_smallestDegreeLast;$ncolors_dinMaxDegree;$ncolors_dSatur;$ncolors_connectedMinDegree,$ncolors_connectedMaxDegree\n")
		end
	end

	return bestcor,bestcolors,bestorder

end


##################################
#Vertex orderings
################################
function remove!(a::Array{Int}, item::Int)
    deleteat!(a, findall(x->x==item, a))
end


function max_aux(elems::Array{Int}, deg::Array{Int})
	v = elems[1]
	max = deg[v]
	for k in elems
		if (deg[k] > max) || (deg[k] == max && k<v)
			v = k
			max = deg[k]
		end
	end
	return v
end

function min_aux(elems::Array{Int}, deg::Array{Int})
	v = elems[1]
	min = deg[v]
	for k in elems
		if (deg[k] < min) || (deg[k] == min && k<v)
			v = k
			min = deg[k]
		end
	end
	return v
end

function sat_aux(elems::Array{Int}, deg::Array{Int}, sat::Array{Int})
	v = elems[1]
	maxsat = sat[v]
	maxdeg = deg[v]
	for k in elems
		if (sat[k] > maxsat) || ((sat[k] == maxsat) && (deg[k] > maxdeg)) || ((sat[k] == maxsat) && (deg[k] == maxdeg) && (k<v))
			v = k
			maxdeg = deg[k]
			maxsat = sat[k]
		end
	end
	return v
end

function maxDegree(inst::InstanceData)
	NV = inst.NV
	G = inst.inputG

	elems = zeros(Int, NV, 2)
	for i = 1:NV
		elems[i,2] = i
		elems[i,1] = length(G[i])
	end

	sort!(elems, dims=1, rev=true)
	order = zeros(Int, NV)

	for i = 1:NV
		order[i] = elems[i,2]
	end

	return order
end

function connectedMaxDegree(inst::InstanceData)
	NV = inst.NV
	G = inst.inputG

	NV = inst.NV
	G = inst.inputG

	visited = zeros(Int,NV)
	order = zeros(Int, NV)

	Q = Vector{Int}()
	degs = [length(G[i]) for i = 1:NV]

	coloredvertices = 0
	while coloredvertices < inst.NV
		v = findmax_aux(degs,visited)
		visited[v] = 1
		coloredvertices += 1
		order[coloredvertices] = v
		for u in G[v]
			if visited[u] == 0
				push!(Q,u)
				visited[u] = 1
			end
		end
		
		while !isempty(Q)
			v = max_aux(Q,degs)
			coloredvertices += 1
			order[coloredvertices] = v
			for u in G[v]
				degs[u] = degs[u] - 1
				if visited[u] == 0
					push!(Q,u)
					visited[u] = 1
				end
			end
			remove!(Q, v)
		end
	end

	return order
end

function dinMaxDegree(inst::InstanceData)
	NV = inst.NV
	G = inst.inputG

	elems = [i for i = 1:NV]
	degs = [length(G[i]) for i = 1:NV]

	order = zeros(Int, NV)
	i =1
	while !isempty(elems)
		v = max_aux(elems,degs)
		order[i] = v
		i = i+1
		remove!(elems, v)
		for u in G[v]
			degs[u] = degs[u] - 1
		end
	end

	return order
end

function minDegree(inst::InstanceData)
	NV = inst.NV
	G = inst.inputG

	elems = [i for i = 1:NV]
	degs = [length(G[i]) for i = 1:NV]

	order = zeros(Int, NV)
	i =1
	while !isempty(elems)
		v = min_aux(elems,degs)
		order[i] = v
		i = i+1
		remove!(elems, v)
		for u in G[v]
			degs[u] = degs[u] - 1
		end
	end

	return order
end

function connectedMinDegree(inst::InstanceData)
	NV = inst.NV
	G = inst.inputG

	visited = zeros(Int,NV)
	order = zeros(Int, NV)

	Q = Vector{Int}()
	degs = [length(G[i]) for i = 1:NV]


	coloredvertices = 0

	while coloredvertices < inst.NV
		v = findmin_aux(degs,visited)

		visited[v] = 1
		coloredvertices += 1
		order[coloredvertices] = v
		for u in G[v]
			if visited[u] == 0
				push!(Q,u)
				visited[u] = 1
			end
		end

		while !isempty(Q)
			v = min_aux(Q,degs)
			coloredvertices += 1
			order[coloredvertices] = v
			for u in G[v]
				degs[u] = degs[u] - 1
				if visited[u] == 0
					push!(Q,u)
					visited[u] = 1
				end
			end
			remove!(Q, v)
		end
	end

	return order
end

function dsatur(inst::InstanceData)
	NV = inst.NV
	G = inst.inputG

	elems = [i for i = 1:NV]
	degs = [length(G[i]) for i = 1:NV]
	sats = [0 for i = 1:NV]
	cor  = 0
	colors = [0 for i = 1:NV]

	order = zeros(Int, NV)
	i = 1
	while !isempty(elems)
		v = sat_aux(elems,degs,sats)
		order[i] = v
		i = i+1
		colors[v] = getMinColor(G[v], colors, v, cor)
		if colors[v] > cor
			cor = colors[v]
		end
		remove!(elems, v)
		for u in G[v]
			degs[u] = degs[u] - 1

			addng = 1
			for z in G[u]
				if colors[v] == colors[z]
					addng = 0
				end
			end
			sats[u] = sats[u] + addng
		end
	end

	return order
end

function getMinColor(neig::Vector{Int},color::Array{Int}, v::Int, cor::Int)

	usedcolors = zeros(Int,cor+1)

	for u in neig
		if color[u] > 0
			usedcolors[color[u]] += 1
		end
	end

	mincolor = -1
	for col in 1:length(usedcolors)
		if usedcolors[col] == 0
			mincolor = col
			break
		end
	end

	return mincolor

end

function getcoloring(G::Array{Vector{Int}}, order::Array{Int}, n::Int)
	cor = 0
	color = [-1 for i = 1:n]

	for i = 1:n
		v = order[i]
		color[v] = getMinColor(G[v], color, v, cor)
		if color[v] > cor
			cor = color[v]
		end
	end
	return (cor, color)
end

########################################################
########################################################
########################################################

function getmincolor(inst::InstanceData, params::ParameterData,v::Int,sol::Array{Int})

	usedcolors = zeros(Int,inst.maxdegree+1)

	for u in inst.inputG[v]
		if sol[u] > 0
			usedcolors[sol[u]] += 1
		end
	end

	mincolor = -1
	for col in 1:inst.NV
		if usedcolors[col] == 0
			mincolor = col
			break
		end
	end

	return mincolor
end

function findmin_aux(degs::Array{Int},visited::Array{Int})
	mindegree = length(visited)+1
	v = -1
	for u in 1:length(visited)
		if visited[u] == 0 && degs[u] < mindegree
			mindegree = degs[u]
			v = u
		end
	end

	return v
end #findmin_aux

function findmax_aux(degs::Array{Int},visited::Array{Int})
	maxdegree = -1
	v = -1
	for u in 1:length(visited)
		if visited[u] == 0 && degs[u] > maxdegree
			maxdegree = degs[u]
			v = u
		end
	end

	return v
end #findmin_aux

end #module