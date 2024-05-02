module Formulations

using JuMP
using Gurobi
using Data
using Solution
using Parameters
using DataStructures

using Printf

#Auxiliar Function to get a Gurobi Model
function defineModel(inst::InstanceData, params::ParameterData)

	model = Model(Gurobi.Optimizer)
	
	#Setting General Parameters
	set_optimizer_attribute(model, "TimeLimit", params.maxtime)
	set_optimizer_attribute(model, "MIPGap", params.tolgap)
	set_optimizer_attribute(model, "Threads", 1)
	
	if params.disablesolver == 1 #Disable gurobi cuts and presolve
		set_optimizer_attribute(model, "Cuts", params.cuts)
		#set_optimizer_attribute(model, "Presolve", params.maxtime)
		#set_optimizer_attribute(model, "Heuristics", params.maxtime)
	end

	if params.modelConfig == 0
		set_optimizer_attribute(model, "Threads", 6)
	elseif params.modelConfig == 1
		set_optimizer_attribute(model, "Threads", 1)
		set_optimizer_attribute(model, "Method", 2)
	elseif params.modelConfig == 2
		set_optimizer_attribute(model, "Threads", 6)
		set_optimizer_attribute(model, "Method", 2)
	end
	
	return model

end

function stdFormulation(inst::InstanceData, params::ParameterData,colors::Array{Int}, order::Array{Int}, 
	zeta::Int, psi)
	if params.problem == "gcol"
		stdFormulationGrundyColoring(inst::InstanceData, params::ParameterData, colors, zeta, psi)
	end

end

function repFormulation(inst::InstanceData, params::ParameterData, colors::Array{Int},order::Array{Int}, 
	zeta::Int, psi)
	if params.problem == "gcol"
		repFormulationGrundyColoring(inst::InstanceData, params::ParameterData, colors, zeta, psi)
	end

end

function stdFormulationGrundyColoring(inst::InstanceData, params::ParameterData, colors::Array{Int}, zeta, psi)
	println("Running Formulations.stdFormulation")
	model = defineModel(inst, params)
	limitColorPerVertex , verticesCanBeColoredK = generateAuxiliaryStdFormulationStructs(inst, zeta, psi)

	maxdegree = inst.maxdegree
	println("maxdegree = ",maxdegree)
	UB = maxdegree+1
	UB = min(UB, zeta, maximum(psi))

	### Defining variables ###
	@variable(model, w[k=1:UB], Bin)
	@variable(model, x[v=1:inst.NV,k=limitColorPerVertex[v]], Bin)

	### Objective function ###
	@objective(model, Max, sum(w[k] for k=1:UB))

	@constraint(model,
		independentset[edge=1:size(inst.edgeslist,1), k=intersect(limitColorPerVertex[(inst.edgeslist[edge,1])], limitColorPerVertex[(inst.edgeslist[edge,2])])],
   		x[(inst.edgeslist[edge,1]),k] + x[(inst.edgeslist[edge,2]),k] <= 1)
	
	@constraint(model,
		iscovered[v=1:inst.NV],
   		sum(x[v,k] for k in limitColorPerVertex[v]) == 1)

	@constraint(model,
		onlyifused[k=1:UB],
		w[k] <= sum(x[v,k] for v=1:inst.NV if k in limitColorPerVertex[v]))

	
	@constraint(model,
		adjcolors[v=1:inst.NV,k=1:(limitColorPerVertex[v][length(limitColorPerVertex[v])]-1),kp=k+1:limitColorPerVertex[v][length(limitColorPerVertex[v])]],
		x[v,kp] <= sum(x[u,k] for u in intersect(inst.inputG[v], verticesCanBeColoredK[k]))
		)
	
	@constraint(model,
		removesymmetry[k=1:UB-1],
		w[k] >= w[k+1])

	if params.relaxation == 1
		for i in 1:UB
			setcategory(y[j,t],:Cont)
		end
	end
	
	if params.initsol == 1
		kcolor = [0 for k in 1:UB]
		for v in 1:inst.NV, k in limitColorPerVertex[v]
			if colors[v] == k
				kcolor[k] = 1
				set_start_value(x[v,k],1)
			else
				set_start_value(x[v,k],0)
			end
		end

		for k in 1:UB
			if kcolor[k] == 1
				set_start_value(w[k],1)
			else
				set_start_value(w[k],0)
			end
		end	
	end


	t1 = time_ns()
	status = optimize!(model)
	t2 = time_ns()
	elapsedtime = (t2-t1)/1.0e9

	nsolutions = result_count(model)
	
	println("Status = ",status)
	if nsolutions > 0

		bestsol = objective_value(model)
		bestbound = objective_bound(model)
		solvertime = solve_time(model)
		solvernodes = node_count(model)
		
		println("Resultscount = ",nsolutions)
		
		println("bestsol = ", bestsol)
		println("Elapsed = ", elapsedtime)
		println("Solver time = ",solvertime)

		soldata = getSolutionFromSTD(x, w, UB, inst.NV, limitColorPerVertex)
		print(soldata)
		validate = validateGrundySolution(soldata, inst.inputG)

		open("modelLP.txt", "w") do f
    		print(f, model)
		end

		open(params.outputfile,"a") do f
			write(f,";$(params.problem);$(params.gamma);$(params.form);$bestsol;$bestbound;$elapsedtime;$solvertime;$solvernodes;$(params.method);$(params.relaxation);$(params.disablesolver);$(validate);$(status);$(params.modelConfig)\n")
		end
		
	else
		bestsol = "NA"
		bestbound = objective_bound(model)
		solvertime = solve_time(model)
		solvernodes = node_count(model)
		open(params.outputfile,"a") do f
			write(f,";$(params.problem);$(params.gamma);$(params.form);$bestsol;$bestbound;$elapsedtime;$solvertime;$solvernodes;$(params.method);$(params.relaxation);$(params.disablesolver);$(params.modelConfig)\n")
		end
	end

end #function stdFormulation()


#Representativa Formulation
function repFormulationGrundyColoring(inst::InstanceData, params::ParameterData, colors::Array{Int}, zeta, psi)
	println("Running Formulations.repFormulation")
	model = defineModel(inst, params)

	maxdegree =  inst.maxdegree
	println("maxdegree = ",maxdegree)
	UB = min(zeta, maximum(psi))
	
	#anti neigbourhoods
	antiG = Array{Vector{Int}}(undef,inst.NV)
	for v in 1:inst.NV
		antiG[v] = []
	end
	for u in 1:inst.NV, v in u+1:inst.NV
		if inst.incmatrix[u,v] == 0
			push!(antiG[u],v)
			push!(antiG[v],u)
		end
	end

	antiGcol = Array{Vector{Int}}(undef,inst.NV)
	for v in 1:inst.NV
		antiGcol[v] = copy(antiG[v])
		push!(antiGcol[v], v)
	end


	### Defining variables ###
	@variable(model, x[u=1:inst.NV,v=u:inst.NV;u==v || u in antiG[v]], Bin)
	@variable(model, y[u=1:inst.NV,v=1:inst.NV;u!=v], Bin)
	@variable(model, 0 <= phi[v=1:inst.NV] <= min(zeta, maximum(psi)))

	### Objective function ###
	@objective(model, Max,
		sum(x[v,v] for v=1:inst.NV))

	#To avoid compute new lists, using a binary matrix
	isvisited = falses(inst.NV,inst.NV)

	### Adding restrictions ###
	#satble set restrictions
	@constraint(model,
		independentset[u=1:inst.NV,v in antiG[u], w in antiG[u]; u<v && v < w && w in inst.inputG[v]],
		x[u,v] + x[u,w] <= x[u,u])

	for u=1:inst.NV,v in antiG[u], w in antiG[u]; u<v && v < w && w in inst.inputG[v]
		isvisited[u,v] = true
		isvisited[u,w] = true
	end

   	@constraint(model,
		independentsetNew[u=1:inst.NV,v in antiG[u]; u<v],
   		x[u,v] <= x[u,u])

	#All vertex has a color
	@constraint(model,
		iscovered[v=1:inst.NV],
   		sum(x[u,v] for u in antiG[v] if u<v ) + x[v,v] == 1)

	#Grundy color restriction
	@constraint(model,
		grundypropuequalsv[u=1:inst.NV,v in antiGcol[u],p=1:inst.NV;u <=v && u!=p],
		x[u,v] <= sum(x[p,w] for w in inst.inputG[v] if w in antiGcol[p] && p <= w) + 1 - y[p,u])

	#ordering
	@constraint(model,
		enforceprec[u in 1:inst.NV, v in 1:inst.NV; u!=v],
		y[v,u] + y[u,v] >=  x[u,u] + x[v,v] - 1)

	@constraint(model,
		precifrepreA[u in 1:inst.NV, v in 1:inst.NV; u!=v],
		y[u,v] + y[v,u] <=  x[u,u])

	@constraint(model, mtz[u=1:inst.NV, v=1:inst.NV; u!=v],
	phi[u] - phi[v] + 1 <= min(zeta, maximum(psi)) * (1 - y[u,v]))

	#relaxing
	if params.relaxation == 1
		for i in 1:UB
			setcategory(y[j,t],:Cont)
		end
	end

	#including initial solution
	if params.initsol == 1

		maxcol = maximum(colors)
		vwithcolor = Array{Vector{Int}}(undef,maxcol)
		for k in 1:maxcol
			vwithcolor[k] = []
		end
		for v in 1:inst.NV
			push!(vwithcolor[colors[v]],v)
		end

		for k in 1:maxcol
			minindex = minimum(vwithcolor[k])
			for v in vwithcolor[k]
				set_start_value(x[minindex,v],1)
			end
		end

		for k1 in 1:maxcol, k2 in 1:maxcol
			minindex1 = minimum(vwithcolor[k1])
			minindex2 = minimum(vwithcolor[k2])
			if k1 < k2
				set_start_value(y[minindex1, minindex2], 1)
			elseif k1 > k2
				set_start_value(y[minindex1, minindex2], 0)
			end
		end
	end

	t1 = time_ns()
	status = optimize!(model)
	t2 = time_ns()
	elapsedtime = (t2-t1)/1.0e9

	status = termination_status(model)
	nsolutions = result_count(model)

	println("Status = ",status)
	if nsolutions > 0

		bestsol = objective_value(model)
		bestbound = objective_bound(model)
		solvertime = solve_time(model)
		solvernodes = node_count(model)
		
		println("Resultscount = ",nsolutions)
		
		println("bestsol = ", bestsol)
		println("Elapsed = ", elapsedtime)
		println("Solver time = ",solvertime)

		soldata = getSolutionFromREP(x, y, phi, antiG, inst.NV)
		print(soldata)
		validate = validateGrundySolution(soldata, inst.inputG)

		open(params.outputfile,"a") do f
			write(f,";$(params.problem);$(params.gamma);$(params.form);$bestsol;$bestbound;$elapsedtime;$solvertime;$solvernodes;$(params.method);$(params.relaxation);$(params.disablesolver);$(validate);$(status);$(params.modelConfig)\n")
		end
		
	else
		bestsol = "NA"
		bestbound = objective_bound(model)
		solvertime = solve_time(model)
		solvernodes = node_count(model)
		open(params.outputfile,"a") do f
			write(f,";$(params.problem);$(params.gamma);$(params.form);$bestsol;$bestbound;$elapsedtime;$solvertime;$solvernodes;$(params.method);$(params.relaxation);$(params.disablesolver);$(params.modelConfig)\n")
		end
	end
	

end #repFormulationGrundyColoring


function generateAuxiliaryStdFormulationStructs(inst, zeta, psi)

	limitColorPerVertex = Array{Vector{Int}}(undef,inst.NV)
	verticesCanBeColoredK = Array{Vector{Int}}(undef,inst.maxdegree+1)
	
	for v in 1:inst.NV
		limitColorPerVertex[v] = [i for i = 1:min(length(inst.inputG[v])+1,zeta, psi[v])]
	end

	for k in 1:(inst.maxdegree+1)
		verticesCanBeColoredK[k] = []
		for v in 1:inst.NV
			if k in limitColorPerVertex[v]
				push!(verticesCanBeColoredK[k], v)
			end
		end
	end
	
	return limitColorPerVertex, verticesCanBeColoredK
end

end
