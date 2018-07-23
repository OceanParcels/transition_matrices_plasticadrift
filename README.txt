Createrawmatrix.py is used for creating the crossing matrices, frac_coast_model.py for creating the beaching model used, createtransitmatrix.py for converting the crossing matrices to transit matrices and global_fwd_csv(_Netherlands).py for using these transit matrices for computing dispersion in time.

The methodology and use of these files is explained in the Method section of the accompanying paper. 
Adj_coast_model.py is not used to produce results in this paper, but is an example of a different beaching model as explained in the Discussion section of the paper.

All initial values are set in createrawmatrix.py and consequently called in the other files, so e.g. changing the model resolution should be done by setting it in createrawmatrix.py.

The input data used for these files and computed results discussed in the paper can be found here: https://doi.org/10.5281/zenodo.1265457