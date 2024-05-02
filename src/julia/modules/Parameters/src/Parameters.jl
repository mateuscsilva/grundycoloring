module Parameters
#
struct ParameterData
    instName::String
    problem::String
    method::String ### exact
    form::String
    solver::String
    maxtime::Int
    tolgap::Float64
    printsol::Int
    disablesolver::Int
    relaxation::Int
    outputfile::String
    initsol::Int
    separationall::Int
    maxmaxcliques
    germansep::Int
    maxnodes::Int
    disablecliques::Int
    germancliques::Int
    gamma::Float64
    heurseed::Int
    maxsolutions::Int
    upperboundDepth::Int
    cuts::Int
    modelConfig::Int
    outputBound::String
end




export ParameterData, ParameterDataBranchAndCut, ParameterDataHFLeuristic,readInputParameters


function readInputParameters(ARGS)

    println("Running Parameters.readInputParameters")

    ### Set standard values for the parameters ###
    instName="../../instances/dimacs_col/myciel5.col"
    problem="gcol"
    method="upperboundnew"
    form="rep"
    solver = "Gurobi"
    maxtime = 3600
    tolgap = 0.000001
    printsol = 0
    disablesolver = 0
    relaxation = 0
    outputfile = "saidaTeste.txt"
    initsol = 0
    separationall = 0
    maxmaxcliques = 500
    germansep = 1
    maxnodes = 10000
    disablecliques = 0
    germancliques = 0
    gamma = 0.95
    heurseed = 1000
    maxsolutions = 1000
    upperboundDepth = 5
    cuts = 2
    modelConfig = 3
    outputBound = "upperbound.txt"

    ### Read the parameters and set correct values whenever provided ###
    for param in 1:length(ARGS)
        if ARGS[param] == "--method"
            method = ARGS[param+1]
            param += 1
        elseif ARGS[param] == "--inst"
            instName = ARGS[param+1]
            param += 1
        elseif ARGS[param] == "--form"
            form = ARGS[param+1]
            param += 1
        elseif ARGS[param] == "--solver"
            solver = ARGS[param+1]
            param += 1
        elseif ARGS[param] == "--maxtime"
            maxtime = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--tolgap"
            tolgap = parse(Float64,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--printsol"
            printsol = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--disablesolver"
            disablesolver = parse(Float64,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--relaxation"
            relaxation = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--outputfile"
            outputfile = ARGS[param+1]
            param += 1
        elseif ARGS[param] == "--problem"
            problem = ARGS[param+1]
            param += 1
        elseif ARGS[param] == "--initsol"
            initsol = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--separationall"
            separationall = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--maxmaxcliques"
            maxmaxcliques = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--germansep"
            germansep = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--maxnodes"
            maxnodes = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--disablecliques"
            disablecliques = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--germancliques"
            germancliques = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--gamma"
            gamma = parse(Float64,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--heurseed"
            heurseed = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--maxsolutions"
            maxsolutions = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--upperboundDepth"
            upperboundDepth = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--modelConfig"
            modelConfig = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--outputbound"
            outputBound = ARGS[param+1]
            param += 1
        end
    end

    params = ParameterData(instName,problem,method,form,solver,maxtime,tolgap,printsol,disablesolver,relaxation,outputfile,initsol,separationall,maxmaxcliques,germansep,maxnodes,disablecliques,germancliques,gamma,heurseed,maxsolutions, upperboundDepth, cuts, modelConfig, outputBound)

    return params

end ### end readInputParameters


end ### end module
