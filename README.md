# FAIR Simulations

This repository contains the code to generate events for the FAIR 30 GeV experiment.
The code is written in Julia and uses `ThreeBodyDecay.jl` for dynamic model of the decay, as well as `FourVectors.jl` for the handling four-vectors.

## Content

- [`notebooks/generator.jl`](notebooks/generator.jl): Pluto notebook to generate events.
- [`notebooks/XtoppJpsi.jl`](notebooks/XtoppJpsi.jl): Dalitz plot for pp → p p Jψ

## Running

1. Clone the repository

```
bash> git clone https://github.com/mmikhasenko/FairSim.git
```

2. Start julia and 

```
bash> julia
```

3. install the dependencies

```julia
] activate .
] instantiate
```

4. Run the code with Pluto

```julia
julia> using Pluto
julia> Pluto.run()
```

5. Open the notebook `notebooks/generator.jl` and run the code.


## Contributing

If you want to contribute to this project, please open an issue or a pull request.

## License

This code is released under the MIT License.