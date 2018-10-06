# DAMICDiffusion

Some of the intial python code use for simulating the diffusion of the charged particles in silicon by incorporating Coulomb repulsion between charges within the charge clouds.

## Use

The useful functions are written in python and located in the diffusion_model_numeric_solutions.py file. There are several different models available (diffusion only--called Groom after the author on the paper, diffusion and coulomb, and piecewise simulation). Each model has a DiffusionFunc(), which describes the partial differential equation being solved, and a DiffusionEquation(), which sets all the intial conditions and solves the PDE using pythons odeint() function.

Each function returns the radial probability charge density of the cloud as a function of time (NxM matrix where N is the number of time steps and M is the number of spatial steps), the radius axis (Mx1 array), and the time axis (Nx1) array.

After the diffusion models are several analysis functions that were used to investigate different aspects of the models and make plots to view them. You can use these as guides of how the DiffusionEquation() are used. Do note that it is possible the functional structure of the DiffusionEquation() for the different models may have changed, so not all analysis functions may be valid.
