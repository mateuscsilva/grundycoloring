module UpperBounds

using Data
using Parameters

export UBData, zeta,  gNumber, gNumberDepth, degSequence

struct UBData
  N::Int #number of vertices
  degreelist::Vector{Tuple{Int,Int}} #list of d^L degree
  degreelistAUX::Vector{Int} #list of d^L degree
  
  function UBData(G::InstanceData)
  	N = G.NV
  	degreelistAUX = [0 for i in 1:N]
  	
  	degreelist = Vector{Tuple{Int,Int}}()
	for v in 1:N
		push!(degreelist,(v,length(G.inputG[v])))
	end
	sort!(degreelist,by = last,rev = false)
	
	return new(N,degreelist, degreelistAUX)
  end
  
end

function calculateUpperBounds(inst::InstanceData,params::ParameterData)
	print("Calculating upper bounds")

	bound_trivial = inst.maxdegree + 1
	bound_delta_two = deltaTwo(inst) + 1
	bound_zeta = newZeta(inst)
	bound_degSeq = degSequence(inst, params)

	open(params.outputBound, "a") do f
		write(f,"$(params.instName);$(bound_trivial);$(bound_zeta);$(bound_delta_two);$(maximum(bound_degSeq))\n")
	end

end

function deltaTwo(G::InstanceData)
	deg = []
	for v in 1:G.NV
		aux = []
		for u in G.inputG[v]
			if length(G.inputG[u]) <= length(G.inputG[v])
				push!(aux, length(G.inputG[u]))
			end
		end
		if(!isempty(aux))
			push!(deg, maximum(aux))
		else 
			push!(deg, 0)
		end
	end
	return maximum(deg)
end

#check the neighboorhood of G for the proper value of k
function checkNeg(G::InstanceData, data::UBData, v::Int)
	k = 0
	for i in 1:G.NV
		if (i in G.inputG[v]) && (data.degreelist[data.degreelist[i][1]][2] > k)
			k = k + 1
		end
	end
	data.degreelistAUX[v] = k
end

function updateData(data::UBData)
	for i in 1:data.N
		data.degreelist[i] = (i, data.degreelistAUX[i])
	end
	sort!(data.degreelist,by = last,rev = false)
end

function zeta(G::InstanceData, l::Int)
	S = []
	k = 0
	data = UBData(G)
	
	for j in 1:l
		for v in 1:G.NV
			checkNeg(G,data,v)
		end
		updateData(data)
	end
	
	for i in 1:G.NV
		if data.degreelist[data.degreelist[i][1]][2] > k
			push!(S, data.degreelist[i][1])
			k = k + 1
		end
	end
	
	return (S,k)
end


function newZeta(inst::InstanceData)
	
	G = deepcopy(inst.inputG)
	residualdegrees = zeros(Int,inst.NV)

	degrees = Vector{Tuple{Int,Int}}()
	for v in 1:inst.NV
		push!(degrees,(length(G[v]),v))
	end
	
	for i in 1:inst.NV
		maxind = argmax(degrees)
		maxv = degrees[maxind][2]
		residualdegrees[i] = degrees[maxind][1]

		deleteat!(degrees,maxind)
		for j in 1:length(degrees)
			(value,key) = degrees[j]
			if key in G[maxv]
				degrees[j] = (degrees[j][1] - 1,degrees[j][2])
			end
		end


	end

	k = inst.maxdegree + 1
	for i in 1:length(residualdegrees)
		if residualdegrees[i]+ i < k
			k = residualdegrees[i]+ i
		end
	end

	return k
end


function gNumber(inst::InstanceData,params::ParameterData)
	NV = inst.NV
	NE = inst.NE

	degree = [length(inst.inputG[v]) for v = 1:NV]
	sort!(degree)
	gnum = 0

	for i = 1:NV
		if degree[i] >= gnum
			gnum = gnum + 1
		end 
	end
	
	return gnum
end

function calculateMaximumGnumOfVertex(inst::InstanceData, vertex::Int, actualLevel::Int, endLevel::Int, path::Array{Int})
	gnum = 1
	actualDegree = 1
	degree = [(length(inst.inputG[v]), v) for v = inst.inputG[vertex] if v ∉ path]
	sort!(degree)

	for i in 1:length(degree)
		aux_degree = degree[i][1] + 1
		if(actualLevel < endLevel)
			push!(path, vertex)
			aux_degree = calculateMaximumGnumOfVertex(inst, degree[i][2], actualLevel+1, endLevel, path)
			deleteat!(path, length(path))
		end	
		
		if aux_degree >= actualDegree
			actualDegree = actualDegree + 1
			gnum = gnum + 1
		end	
	end	

	return gnum
end

function gNumberDepth(inst::InstanceData,params::ParameterData)
	NV = inst.NV
	NE = inst.NE

	degree = [(length(inst.inputG[v]), v) for v = 1:NV]
	gnum = 0
	sort!(degree, rev=true)

	gnum = 1
	for i in 1:NV
		vertexGnum = calculateMaximumGnumOfVertex(inst, degree[i][2], 0, params.upperboundDepth, Int[])
		if vertexGnum > gnum
			gnum = vertexGnum
		end
	end
	
	return gnum
	
end



function degSequence(inst::InstanceData,params::ParameterData)
	H = zeros(Int,inst.NV, inst.maxdegree+1)

	for v in 1:inst.NV
		H[v,1] = 1
	end

	for k in 2:inst.maxdegree+1, v in 1:inst.NV
		H[v,k] = calculateL(inst,v,k,H)
	end
	return H[:,inst.maxdegree+1]
end

function calculateL(inst::InstanceData,v::Int,k::Int,H::Matrix{Int64})
	vec = zeros(Int,length(inst.inputG[v])+1)
	i=1

	for u in inst.inputG[v]
		vec[i] = H[u,k-1]
		i += 1
	end

	vec[i] = 0
	sort!(vec)
	M = zeros(Int,length(vec))

	for j in 2:length(vec)
		if vec[j] <= M[j-1] 
			M[j] = M[j-1]
		else
			M[j] = M[j-1] + 1
		end
	end

	return (M[length(M)] + 1)

end

end
