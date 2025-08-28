# Gradient Descent for reducing Hopping Terms in a MoS2 Tight Binding model
The folder "graddesc_mos2" contains the main scripts of the project, where Nesterov Gradient Descent algorithms are implemented to reduce the complexity of the MoS2 tight binding model, while keeping a good estimate of the bandstructure.  

The following variations of the algorithm are implemented:

- **mos2_graddesc.py**  
The standard NGD algorithm is implemented here for reducing hoppings. The algorithm optimizes for both first and second degree neighbours at the same time.  

- **mos2_multilevel_graddesc.py**  
In this variation of the algorithm, only first OR second degree neighbours are optimized for at a certain time span. Variations are possible where there are e.g. 200 iterations of optimization for first degree neighbours followed by 300 iterations of optimization for second degree neighbours..

- **mos2_bandgab_graddesc.py**  
Here the algorithm optimizes both degrees of order at the same time agein. However the error metric takes the bandstructure-error near the bandgap into account at a much higher weight, leading to clearer bandstructures near the bandgap.  

## Dependencies
All three scripts are based on the class "mcell" and functions in "resources/mos2class.py".  
Make sure that you have followed all steps in the [installation description](../README.md) of the project correctly.