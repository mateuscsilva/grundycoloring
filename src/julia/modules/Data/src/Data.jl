module Data

using EzXML, GraphIO, LightGraphs

struct InstanceData
  NV::Int #number of vertices
  NE::Int #number of edges
  edgeslist::Matrix{Int}
  incmatrix::Matrix{Int}
  inputG::Array{Vector{Int}}
  inputGedges::Array{Vector{Int}}
  incmatrixedges::Matrix{Int}
  maxdegree::Int#maximum degree
end

export InstanceData, readDataMtx, maxDegreeOrder


#change the graph for the first vertices being the ones with largest degree
function maxDegreeOrder(inst::InstanceData)
 	NV = inst.NV
	NE = inst.NE
	maxdegree = inst.maxdegree

	auxN = [i for i = 1:NV ]
	f = [i for i = 1:NV ]
	degree = zeros(Int,NV)

	for v in 1:NV
		degree[v] = length(inst.inputG[v])
	end
	
	#print(degree)

	sorteddeg = Vector{Tuple{Int,Int}}()
	for v in 1:NV
		push!(sorteddeg,(v,degree[v]))
	end
	sort!(sorteddeg,by = last,rev = true)
	
	#print(sorteddeg)

	for v in 1:NV
		degree[v] = sorteddeg[v][2]
		auxN[v] = sorteddeg[v][1]
	end


	for i in 1:NV
		f[auxN[i]] = i
	end

	edgeslist = zeros(Int,NE,2)

	for edge in 1:NE
		edgeslist[edge, 1] = f[inst.edgeslist[edge,1]]
		edgeslist[edge, 2] = f[inst.edgeslist[edge,2]]
	end

	inputGedges = Array{Vector{Int}}(undef,NV)
	inputG = Array{Vector{Int}}(undef,NV)
	for v in 1:NV
		inputG[v] = []
		inputGedges[v] = []
	end

	for edge in 1:NE
		u = edgeslist[edge,1]
		v = edgeslist[edge,2]
		push!(inputG[u],v)
		push!(inputG[v],u)
		push!(inputGedges[u],edge)
		push!(inputGedges[v],edge)
	end


	

	incmatrix = zeros(Int,NV,NV)
	incmatrixedges = zeros(Int,NV,NV)
	for edge in 1:NE
		incmatrix[edgeslist[edge,1],edgeslist[edge,2]] = 1
		incmatrix[edgeslist[edge,2],edgeslist[edge,1]] = 1

		incmatrixedges[edgeslist[edge,1],edgeslist[edge,2]] = edge
		incmatrixedges[edgeslist[edge,2],edgeslist[edge,1]] = edge
	end


	graphdensity = (NE/(NV*(NV-1)/2))
	println("Graph density = $graphdensity")


	instNEW = InstanceData(NV, NE,edgeslist,incmatrix,inputG,inputGedges,incmatrixedges,maxdegree)

	return instNEW

end

function readDataGraphml(instanceFile,params)

	println("Running Data.readData with file $(instanceFile)")


	D = loadgraphs(instanceFile, GraphIO.GraphML.GraphMLFormat())

	NV=1
	NE=1

	G = get(D, "G", get(D,"graph",[]))


	NV = nv(G)
	NE = ne(G)
	edgeslist = zeros(Int,NE,2) #demand NI x NT


	edge = 0


	graphdensity = (NE/(NV*(NV-1)/2))
	println("Graph density = $graphdensity")

	for e in edges(G)
		edge += 1
		edgeslist[edge,1] = src(e)
		edgeslist[edge,2] = dst(e)
	end


	inputGedges = Array{Vector{Int}}(undef,NV)
	inputG = Array{Vector{Int}}(undef,NV)
	for v in 1:NV
		inputG[v] = []
		inputGedges[v] = []
	end
	for edge in 1:NE
		u = edgeslist[edge,1]
		v = edgeslist[edge,2]
		push!(inputG[u],v)
		push!(inputG[v],u)
		push!(inputGedges[u],edge)
		push!(inputGedges[v],edge)
	end


	incmatrix = zeros(Int,NV,NV)
	incmatrixedges = zeros(Int,NV,NV)
	for edge in 1:NE
		incmatrix[edgeslist[edge,1],edgeslist[edge,2]] = 1
		incmatrix[edgeslist[edge,2],edgeslist[edge,1]] = 1

		incmatrixedges[edgeslist[edge,1],edgeslist[edge,2]] = edge
		incmatrixedges[edgeslist[edge,2],edgeslist[edge,1]] = edge
	end


	maxdegree = 0
	for v in 1:NV
		if length(inputG[v]) > maxdegree
			maxdegree = length(inputG[v])
		end
	end

	inst = InstanceData(NV, NE,edgeslist,incmatrix,inputG,inputGedges,incmatrixedges,maxdegree)

	return inst

end





function readDataMtx(instanceFile,params)

	println("Running Data.readData with file $(instanceFile)")
	file = open(instanceFile)

	lines = readlines(file)

	tokens = split(lines[2])
	NV = parse(Int,tokens[1])
	NE = parse(Int,tokens[length(tokens)])


	edgeslist = zeros(Int,NE,2) #demand NI x NT
	edge = 0
	for l in 3:length(lines)
		tokens = split(lines[l])
		edge += 1
		edgeslist[edge,1] = parse(Int,tokens[1])
		edgeslist[edge,2] = parse(Int,tokens[2])
	end

	inputGedges = Array{Vector{Int}}(undef,NV)
	inputG = Array{Vector{Int}}(undef,NV)
	for v in 1:NV
		inputG[v] = []
		inputGedges[v] = []
	end
	for edge in 1:NE
		u = edgeslist[edge,1]
		v = edgeslist[edge,2]
		push!(inputG[u],v)
		push!(inputG[v],u)
		push!(inputGedges[u],edge)
		push!(inputGedges[v],edge)
	end

	incmatrix = zeros(Int,NV,NV)
	incmatrixedges = zeros(Int,NV,NV)
	for edge in 1:NE
		incmatrix[edgeslist[edge,1],edgeslist[edge,2]] = 1
		incmatrix[edgeslist[edge,2],edgeslist[edge,1]] = 1

		incmatrixedges[edgeslist[edge,1],edgeslist[edge,2]] = edge
		incmatrixedges[edgeslist[edge,2],edgeslist[edge,1]] = edge
	end

	graphdensity = (NE/(NV*(NV-1)/2))
	println("Graph density = $graphdensity")
	G = SimpleGraph(NV);
	for edge in 1:NE
		add_edge!(G, edgeslist[edge,1],edgeslist[edge,2]);
	end

	newinstname = string(params.instName,".graphml")

    maxdegree = 0
    for v in 1:NV
        if length(inputG[v]) > maxdegree
            maxdegree = length(inputG[v])
        end
    end

    inst = InstanceData(NV, NE,edgeslist,incmatrix,inputG,inputGedges,incmatrixedges,maxdegree)

	return inst

end


function readDataCol(instanceFile,params)

	println("Running Data.readDataCol with file $(instanceFile)")
	file = open(instanceFile)

	edgestuples = Vector{Tuple{Int,Int}}()
	
	NV=0
	NE=0
	auxedge = 0
	for ln in eachline(file)
		if length(ln) > 0
			if ln[1] == 'p'
				test = split(ln)
				NV = parse(Int,test[3])
				
			elseif ln[1] == 'e'
				auxedge = auxedge + 1
				test = split(ln)
				v1 = parse(Int,test[2])
				v2 = parse(Int,test[3])
				push!(edgestuples,(min(v1,v2),max(v1,v2)))
			end

		end
	end



	edgestuples = unique(edgestuples)

	NE = length(edgestuples)
	edgeslist = zeros(Int,NE,2) #demand NI x NT
	e=0
	for (v1,v2) in edgestuples
		e += 1
		edgeslist[e,1] = v1
		edgeslist[e,2] = v2
	end

	inputGedges = Array{Vector{Int}}(undef,NV)
	inputG = Array{Vector{Int}}(undef,NV)
	for v in 1:NV
		inputG[v] = []
		inputGedges[v] = []
	end
	for edge in 1:NE
		u = edgeslist[edge,1]
		v = edgeslist[edge,2]
		push!(inputG[u],v)
		push!(inputG[v],u)
		push!(inputGedges[u],edge)
		push!(inputGedges[v],edge)
	end

	for v in 1:NV
		inputG[v] = unique(inputG[v])
	end

	incmatrix = zeros(Int,NV,NV)
	incmatrixedges = zeros(Int,NV,NV)
	for edge in 1:NE
		incmatrix[edgeslist[edge,1],edgeslist[edge,2]] = 1
		incmatrix[edgeslist[edge,2],edgeslist[edge,1]] = 1

		incmatrixedges[edgeslist[edge,1],edgeslist[edge,2]] = edge
		incmatrixedges[edgeslist[edge,2],edgeslist[edge,1]] = edge
	end

	correctedNE = sum(length(inputG[v]) for v in 1:NV)/2
	graphdensity = (correctedNE/(NV*(NV-1)/2))
	open(params.outputfile,"a") do f
    write(f,";$(graphdensity)")
	end
	println("Graph density = $graphdensity")

	G = SimpleGraph(NV);
	for edge in 1:NE
		add_edge!(G, edgeslist[edge,1],edgeslist[edge,2]);
	end

	maxdegree = 0
	for v in 1:NV
		if length(inputG[v]) > maxdegree
			maxdegree = length(inputG[v])
		end
	end

	inst = InstanceData(NV, NE,edgeslist,incmatrix,inputG,inputGedges,incmatrixedges,maxdegree)
	return inst

end


function readDataClq(instanceFile,params)

	println("Running Data.readDataClq with file $(instanceFile)")
	file = open(instanceFile)

	edgestuples = Vector{Tuple{Int,Int}}()
	
	NV=0
	NE=0
	auxedge = 0
	for ln in eachline(file)
		if length(ln) > 0
			if ln[1] == 'p'
				test = split(ln)
				NV = parse(Int,test[3])
			elseif ln[1] == 'e'
				auxedge = auxedge + 1
				test = split(ln)
				v1 = parse(Int,test[2])
				v2 = parse(Int,test[3])
				push!(edgestuples,(min(v1,v2),max(v1,v2)))
			end
		end
	end

	edgestuples = unique(edgestuples)
	NE = length(edgestuples)
	edgeslist = zeros(Int,NE,2) #demand NI x NT
	e=0
	for (v1,v2) in edgestuples
		e += 1
		edgeslist[e,1] = v1
		edgeslist[e,2] = v2
	end

	inputGedges = Array{Vector{Int}}(undef,NV)
	inputG = Array{Vector{Int}}(undef,NV)
	for v in 1:NV
		inputG[v] = []
		inputGedges[v] = []
	end
	for edge in 1:NE
		u = edgeslist[edge,1]
		v = edgeslist[edge,2]
		push!(inputG[u],v)
		push!(inputG[v],u)
		push!(inputGedges[u],edge)
		push!(inputGedges[v],edge)
	end

	for v in 1:NV
		inputG[v] = unique(inputG[v])
	end

	incmatrix = zeros(Int,NV,NV)
	incmatrixedges = zeros(Int,NV,NV)
	for edge in 1:NE
		incmatrix[edgeslist[edge,1],edgeslist[edge,2]] = 1
		incmatrix[edgeslist[edge,2],edgeslist[edge,1]] = 1

		incmatrixedges[edgeslist[edge,1],edgeslist[edge,2]] = edge
		incmatrixedges[edgeslist[edge,2],edgeslist[edge,1]] = edge
	end

	correctedNE = sum(length(inputG[v]) for v in 1:NV)/2
	graphdensity = (correctedNE/(NV*(NV-1)/2))
	println("Graph density = $graphdensity")

	G = SimpleGraph(NV);
	for edge in 1:NE
		add_edge!(G, edgeslist[edge,1],edgeslist[edge,2]);
	end

	maxdegree = 0
	for v in 1:NV
		if length(inputG[v]) > maxdegree
			maxdegree = length(inputG[v])
		end
	end

	inst = InstanceData(NV, NE,edgeslist,incmatrix,inputG,inputGedges,incmatrixedges,maxdegree)

	return inst

end

function readDataBcol(instanceFile,params)

	println("Running Data.readDataBcol with file $(instanceFile)")
	file = open(instanceFile)

	edgestuples = Vector{Tuple{Int,Int}}()
	
	iterations = 0
	NV=0
	NE=0
	auxedge = 0
	for ln in eachline(file)
		if length(ln) > 0
			if iterations == 0
				iterations = 1
				NV = parse(Int,ln)
			else
				auxedge = auxedge + 1				
				test = split(ln)
				v1 = parse(Int,test[1]) + 1
				v2 = parse(Int,test[2]) + 1
				push!(edgestuples,(min(v1,v2),max(v1,v2)))
			end

		end
	end

	edgestuples = unique(edgestuples)

	NE = length(edgestuples)
	edgeslist = zeros(Int,NE,2) #demand NI x NT
	e=0
	for (v1,v2) in edgestuples
		e += 1
		edgeslist[e,1] = v1
		edgeslist[e,2] = v2
	end

	inputGedges = Array{Vector{Int}}(undef,NV)
	inputG = Array{Vector{Int}}(undef,NV)
	for v in 1:NV
		inputG[v] = []
		inputGedges[v] = []
	end
	for edge in 1:NE
		u = edgeslist[edge,1]
		v = edgeslist[edge,2]
		push!(inputG[u],v)
		push!(inputG[v],u)
		push!(inputGedges[u],edge)
		push!(inputGedges[v],edge)
	end

	for v in 1:NV
		inputG[v] = unique(inputG[v])
	end

	incmatrix = zeros(Int,NV,NV)
	incmatrixedges = zeros(Int,NV,NV)
	for edge in 1:NE
		incmatrix[edgeslist[edge,1],edgeslist[edge,2]] = 1
		incmatrix[edgeslist[edge,2],edgeslist[edge,1]] = 1

		incmatrixedges[edgeslist[edge,1],edgeslist[edge,2]] = edge
		incmatrixedges[edgeslist[edge,2],edgeslist[edge,1]] = edge
	end

	correctedNE = sum(length(inputG[v]) for v in 1:NV)/2
	graphdensity = (correctedNE/(NV*(NV-1)/2))
	println("Graph density = $graphdensity")

	G = SimpleGraph(NV);
	for edge in 1:NE
		add_edge!(G, edgeslist[edge,1],edgeslist[edge,2]);
	end

	maxdegree = 0
	for v in 1:NV
		if length(inputG[v]) > maxdegree
			maxdegree = length(inputG[v])
		end
	end
	inst = InstanceData(NV, NE,edgeslist,incmatrix,inputG,inputGedges,incmatrixedges,maxdegree)
	return inst

end


end
