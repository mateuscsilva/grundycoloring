# Formulations

## Run standard formulation

- julia grundycoloring.jl --inst path_to_instance --method exact --form std --maxtime 3600

## Run representative formulation

- julia grundycoloring.jl --inst path_to_instance --method exact --form rep --maxtime 3600

### Other parameters

- tolgap
- disablesolver
- modelConfig
- outputfile

## Package versions

- GraphIO v"0.7.0"
- DataStructures v"0.18.16"
- EzXML v"1.2.0"
- LightGraphs v"1.3.5"
- Gurobi v"1.2.1"
- LightGraphsFlows v"0.3.1"
- JuMP v"1.20.0"