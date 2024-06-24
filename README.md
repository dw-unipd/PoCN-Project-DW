# Physics of Complex Networks 2024 project under Professor De Domenico
Student name: David Weingut

Tasks:

|name|number|score|type|
|-|-|-|-|
|Ising model|1|0.4|theory|
|Traffic congestion|16|0.6|theory|
|Public transport in large cities worldwide|41|1.2|data|


## Running the tasks
Each project can be run after installing the julia programming from [the official website](https://julialang.org/downloads/) and instantiating the projects inside the directories by `julia --project=code/chosenproject -e "import Pkg; Pkg.instantiate()`

for Project #01 additionally the following needs to be run inside the activated julia environment
```julia
using PyCall
run(`$(PyCall.python) -m pip install python-igraph`)
```

#41 should be quick to run while #01 and #16 can take several hours depending on hardware supplied to the julia process(es)