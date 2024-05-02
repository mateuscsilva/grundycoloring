module Solution

using JuMP
using Gurobi

struct SolData
  NV::Int #number of vertices
  NC::Int #number of colors
  colorlist::Array{Int}
  vertexorder::Array{Int} #each position i is the place vertex i has in the order
  ordervertex::Array{Int} #each position i is the vertex that is placed in the position i
end

export SolData, validateGrundySolution, validateGrundyConnectedSolution, getSolutionFromSTD, getSolutionFromREP, getSolutionFromSTDCon, getSolutionFromREPCon

function validateGrundySolution(sol::SolData, inputG::Array{Vector{Int}})
	v = 0

	#all vertices were colored?
	for i = 1:sol.NV
		if sol.colorlist[i] < 0
			println("vertex $i has NO COLOR")
			return false
		end
	end

	neg = zeros(Int64, sol.NC)

	#grundy property
	for i = 1:sol.NV
		v = sol.ordervertex[i]
		for k = 1:sol.NC
			neg[k] = 0
		end

		for u in inputG[v]
			if sol.vertexorder[u] < sol.vertexorder[v]
				neg[sol.colorlist[u]] = neg[sol.colorlist[u]] + 1
			end
		end

		for k = 1:(sol.colorlist[v]-1)
			if neg[k] == 0
				println("vertex $v has NO NEIGHBOOR in color $k")
				return false
			end
		end
	end
	println("SOLUTION OK")
	return true
end

function validateGrundyConnectedSolution(sol::SolData, inputG::Array{Vector{Int}})

	#verify connectivity
	list = []
	push!(list,sol.ordervertex[1])
	for i = 2:sol.NV
		v = sol.ordervertex[i]
		if isempty(intersect(list, inputG[v]))
			println("vertex $v has NO NEIGHBOOR COLORED!")
			return false
		end
		push!(list,v)
	end

	return validateGrundySolution(sol, inputG)
end

function getSolutionFromSTD(x, w, m, n, limitColorPerVertex)
	NV = n
	NC = 0
	colorlist = [-1 for i=1:n]
	vertexorder = [-1 for i=1:n]
	ordervertex = [-1 for i=1:n]

	od = 1
	for k=1:m
		if value(w[k]) >= 0.0001
		#if w[k] >= 0.01
			NC=NC+1
			for v = 1:n
				if k in limitColorPerVertex[v] && value(x[v,k])>= 0.01
				#if x[v,k] >= 0.01
					colorlist[v] = NC
					vertexorder[v] = od
					ordervertex[od] = v
					od = od+1
				end
			end
		end
	end

	return SolData(NV, NC, colorlist, vertexorder, ordervertex)
end

function getSolutionFromREP(x, y, phi, antiG, n::Int)
	NV = n
  NC = 0
  colorlist = [-1 for i=1:n]
  vertexorder = [-1 for i=1:n]
  ordervertex = [-1 for i=1:n]

	vcol = [v for v=1:n if value(x[v,v])>= 0.01]
	NC = length(vcol)

	
	for v = 1:NC
		for u = 1:(NC-1)
			if value(y[vcol[u+1],vcol[u]])>= 0.01
				aux = vcol[u]
				vcol[u] = vcol[u+1]
				vcol[u+1] = aux
			end
		end
	end

	
	#for u in 1:NV
	#	if(value(x[u,u])>= 0.01)
	#		println("vertex $u is a REP")
	#	end
	#	for v in antiG[u]
	#		if (u < v) && (value(x[u,v])>= 0.01)
	#			println("vertex $v has REP $u")
	#		end
	#	end
	#end	

	od = 1
	for k = 1:NC
		u = vcol[k]
		colorlist[u] = k
		vertexorder[u] = od
		ordervertex[od] = u
		od = od+1
		for v in antiG[u]
			if (u < v) && (value(x[u,v])>= 0.01)
				colorlist[v] = k
				vertexorder[v] = od
				ordervertex[od] = v
				od = od+1
			end
		end
	end
	return SolData(NV, NC, colorlist, vertexorder, ordervertex)
end

function getSolutionFromSTDCon(z, w, m::Int, n::Int, limitColorPerVertexWithTime)
	NV = n
	NC = 0
	colorlist = [-1 for i=1:n]
	vertexorder = [-1 for i=1:n]
	ordervertex = [-1 for i=1:n]
	
	od = 1
	for t=1:n
		for v=1:n
			for k in limitColorPerVertexWithTime[v][t]
				if value(z[v,t,k]) >= 0.01
					colorlist[v] = k
					vertexorder[v] = t
					ordervertex[t] = v
					od = od+1
				end
			end
		end
	end

	for k=1:m
		if value(w[k]) >= 0.01
			NC = NC + 1
		end
	end

	#for k=1:m
	#	if value(w[k]) >= 0.01
		#if w[k] >= 0.01
	#		NC=NC+1
	#		for v = 1:n
	#			for t = 1:n
	#				if k in limitColorPerVertex[v] && value(z[v,t,k])>= 0.01
	#					colorlist[v] = NC
	#					vertexorder[v] = t
	#					ordervertex[t] = v
	#					od = od+1
	#				end
	#			end
	#		end
	#	end
	#end
	return SolData(NV, NC, colorlist, vertexorder, ordervertex)
end

function getSolutionFromREPCon(x, y, antiG, n::Int)
	NV = n
  NC = 0
  colorlist = [-1 for i=1:n]
  vertexorder = [-1 for i=1:n]
  ordervertex = [-1 for i=1:n]

	vcol = [v for v=1:n,t = 1:n if value(x[v,v,t])>= 0.01]
	NC = length(vcol)

	#for v = 1:NC
	#	for u = 1:(NC-1)
	#		if v != u
	#			if value(y[vcol[u],vcol[u+1]]) >= 0.01
	#				aux = vcol[v]
	#				vcol[v] = vcol[u]
	#				vcol[u] = aux
	#			end
	#	end
	#	end
	#end

	for v = 1:NC
		for u = 1:(NC-1)
			if v != u 
				if value(y[vcol[u+1],vcol[u]])>= 0.01
					aux = vcol[u]
					vcol[u] = vcol[u+1]
					vcol[u+1] = aux
				end
			end
		end
	end

	#vcol =  reverse(vcol)
	#for u in 1:NV
	#	if(value(x[u,u])>= 0.01)
	#		println("vertex $u is a REP")
	#	end
	#	for v in antiG[u]
	#		if (u < v) && (value(x[u,v])>= 0.01)
	#			println("vertex $v has REP $u")
	#		end
	#	end
	#end	

	od = 1
	for t in 1:NV
		for i in 1:NV
			for j in i:NV
				if (i in antiG[j] || i==j)
					if(value(x[i,j,t]) >= 0.0001)
						k = findall(x->x==i,vcol)
						k = k[1]
						colorlist[j] = k
						vertexorder[j] = od
						ordervertex[od] = j
						od = od+1
					end
				end
			end
		end
	end


	#for k = 1:NC
	#	u = vcol[k]
	#	colorlist[u] = k
	#	vertexorder[u] = od
	#	ordervertex[od] = u
	#	od = od+1
	#	for v in antiG[u]
	#		for t = 1:n
	#			if (u < v) && (value(x[u,v,t])>= 0.01)
	#				colorlist[v] = k
	#				vertexorder[v] = od
	#				ordervertex[od] = v
	#				od = od+1
	#			end
	#		end
	#	end
	#end
	return SolData(NV, NC, colorlist, vertexorder, ordervertex)
end

#n = 5
#w = zeros(Int64, n)
#x = zeros(Int64, n, n)

#w[5] = 1
#w[3] = 1
#w[1] = 1

#x[1,1] = 1
#x[2,3] = 1
#x[3,5] = 1
#x[4,1] = 1
#x[5,3] = 1

#soldata = getSolutionFromSTD(x, w, 5, 5)

#println(w)
#println(x)
#println(soldata)

end
