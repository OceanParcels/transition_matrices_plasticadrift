This repository holds code for creating transition matrices, as described in e.g. [Van Sebille 2014](http://dx.doi.org/10.1016/j.jembe.2014.09.002). See the [PDF of the thesis](https://github.com/OceanParcels/transition_matrices_plasticadrift/blob/master/Simulating%20pathways%20and%20beaching%20effects%20of%20plastic%20originating%20from%20the%20Dutch%20coast.pdf) for more details. 

Createrawmatrix.py is used for creating the crossing matrices, frac_coast_model.py for creating the beaching model used, createtransitmatrix.py for converting the crossing matrices to transit matrices and global_fwd_csv(_Netherlands).py for using these transit matrices for computing dispersion in time.

The methodology and use of these files is explained in the Method section of the accompanying paper. 
Adj_coast_model.py is not used to produce results in this thesis, but is an example of a different beaching model as explained in the Discussion section of the [thesis](https://github.com/OceanParcels/transition_matrices_plasticadrift/blob/master/Simulating%20pathways%20and%20beaching%20effects%20of%20plastic%20originating%20from%20the%20Dutch%20coast.pdf).

All initial values are set in createrawmatrix.py and consequently called in the other files, so e.g. changing the model resolution should be done by setting it in createrawmatrix.py.

The input data used for these files and computed results discussed in the paper can be found here: https://doi.org/10.5281/zenodo.1265457
