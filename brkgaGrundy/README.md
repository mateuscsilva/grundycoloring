# BRKGA 

## Compile and Execute BRKGA

1. make
2. ./brkga-grundy --inst path_to_instance --solver grundyNP --seed $seed --time 300 --pe 0.1 --pm 0.15 --rhoe 0.6 --factor 2


### Possible values for pe pm

- pe is the percentage of the elite population.
- pm is the percentage of the mutant population.
- Any value in the range (0,1)
- We indicate that "pe" is in the range 0.1 and 0.3.
- We indicate that "pm" is in the range 0.1 and 0.25.

### Possible values for rhoe

- Indicates the probability of inheritance.
- By definition "rhoe" > 0.5.

### Possible values for time

- Any integer value > 0. 
- Time is in seconds.
- We use the default timeout of 300 seconds.

### Possible values for factor

- Population size is defined as chromosome size times "factor".
- Any integer number > 0.
- The most common options are 1, 2, 3.

### Possible values for seed

- Any integer number