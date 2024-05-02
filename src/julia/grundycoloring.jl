push!(LOAD_PATH, "modules/")
using Pkg

using Gurobi

import Data
import Parameters
import Formulations
import Heuristics
import Enumerative
import Solution
import UpperBounds



# Read the parameters from command line
params = Parameters.readInputParameters(ARGS)

open(params.outputfile,"a") do f
    write(f,"$(params.instName)")
end


# Read instance data
if occursin(".graphml", params.instName)
    inst = Data.readDataGraphml(params.instName,params)
elseif occursin(".mtx", params.instName)
    inst = Data.readDataMtx(params.instName,params)
elseif occursin(".col", params.instName)
    inst = Data.readDataCol(params.instName,params)
elseif occursin(".clq", params.instName)
    inst = Data.readDataClq(params.instName,params)
elseif occursin(".bcol", params.instName)
    inst = Data.readDataBcol(params.instName,params)
end


if params.method == "exact"
    cor,colors,order = Heuristics.coloringHeuristics(inst,params, params.problem == "cgcol")
    zeta = UpperBounds.newZeta(inst)
    deltaTwo = UpperBounds.deltaTwo(inst) + 1
    degSequence = UpperBounds.degSequence(inst, params)
    if params.form == "std"
        Formulations.stdFormulation(inst, params, colors, order, min(zeta,deltaTwo), degSequence)
    elseif params.form == "rep"
        Formulations.repFormulation(inst, params, colors, order, min(zeta,deltaTwo), degSequence)
    end

elseif params.method == "coloringheuristics"
    cor,colors,order = Heuristics.coloringHeuristics(inst,params)
elseif params.method == "upperbounds"
    UpperBounds.calculateUpperBounds(inst, params)    
elseif params.method == "upperbound"
    UpperBounds.gNumber(inst, params)
elseif params.method == "upperboundnew"
    UpperBounds.zeta(inst, params.upperboundDepth)    
end
