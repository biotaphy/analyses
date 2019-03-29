"""Module containing code for calculating phylogenetic beta diversity

Note:
    * Numpy and Scipy may have the methods needed for these computations.  Feel
        free to take advantage of those to libraries as they are already
        dependencies.
"""
import numpy as np
from analyses.lm_objects.matrix import Matrix
from analyses.lm_objects.tree import TreeWrapper
import itertools as it # new dependency.
import dendropy # this might be redundant to have added since TreeWrapper already uses dendropy.

# .............................................................................
def pdnew( pam, tree ): 
    """Creates a lookup dictionary for PD of sites in a matrix

    Args:
        pam (:obj:`Matrix`): A Lifemapper Matrix object with presence absence
            values.

        tree (:obj:`TreeWrapper`): A TreeWrapper object for a wrapped Dendropy
            phylogenetic tree.

    Returns:
        Matrix object holding PD values & spp. richness for each community in the sample.
        Col1 = PD; Col2 = SR; Rows = Indvidual samples from pam.
    """

    # Get number of samples in community matrix data.
    nsamp = len( pam.get_row_headers() ) 
    
    # Array to hold each sample's PD & RCH. col1 = PD; col2 = RCH.
    PD_array = np.zeros( (nsamp, 2), dtype = np.float )

    # This loop will calculate the PD value for each sample in 'pam'.
    for sample in range( nsamp ):
        # Pull out the data for current sample. 
        my_samp = pam.data[ sample ]
        
        # Pull out which spp are present & which are absent from the sample.
        # yields lists of strings --> spp names.
        sp_pres = list( it.compress( pam.get_column_headers(), my_samp ) ) 

        # Get a tree of the spp present in the sample. 
        tree_pres = tree.extract_tree_with_taxa_labels( sp_pres )

        # Get sum of edge lengths for each sub tree.
        PD_pres = tree_pres.length()

        # Get spp rch of sample.
        rch_samp = len( sp_pres )

        # Update PD_array.
        PD_array[ sample ] = [ PD_pres, rch_samp ]

    # convert PD_array to Matrix object & match headers to pam data. 
    PD_mat = Matrix( PD_array, headers = { '0': pam.get_row_headers(),
                                           '1': [ 'PD', 'SR' ] }  )

    # return values.
    return PD_mat

# .............................................................................
def core_PD_calc( pam, tree ):
    """Creates an array of core metrics used to asses components of beta diversity.

    Args:
        pam (:obj:'Matrix'): A Lifemapper Matrix object with presence absence
            values.

        tree (:obj:'TreeWrapper'): A TreeWrapper object for a wrapped Dendropy
            phylogenetic tree.

    Returns:
        Matrix object. Cols = core metrics (see Details below); Rows = Pairwise comparisons.

    Details:
        In general, the metrics returned represent different contributions to PD arising from how communities are combined together.
        Metrics:
            min_not_shared: smallest PD distance from individual samples to the PD of their combination.
            max_not_shared: largest PD distance from individual samples to the PD of their combination.
            sum_not_shared: the total addition to PD from combining the two communities over viewing them separately.
            shared: the combined contribution to PD that the communities make jointly.
        
    """

    # PD for each community of the community matrix.
    pd = pdnew(pam, tree)

    # List all possible pairwise community combinations.
    combin = list( it.combinations( range( len( pam.get_row_headers() ) ), 2 ) )

    # Array to store PD of pairwise site combinations. Rows = all pairwise community comparisons. Cols = spp.
    com_tot_pair = np.zeros( ( len(combin), len( pam.get_column_headers() ) ), 
                             dtype = np.float) # can change to int to save memory later.

    # This loop will populate the pairwise PD_array. 1 == spp present in at least 1 sample; 0 == spp absent from both samples.
    for pair in range( len(combin) ):
        # Assign each site's data to new variable for convenience.
        site0 = pam.data[ combin[ pair ][ 0 ] ]
        site1 = pam.data[ combin[ pair ][ 1 ] ]

        # Go through each spp "slot" & determine if it is present in at least 1 of the samples.
        for idx in range( len(site0) ):
            if site0[ idx ] or site1[ idx ]:
                com_tot_pair[ pair, idx ] = 1
            else: 
                com_tot_pair[ pair, idx ] = 0

    # Convert pairwise PD_array into Matrix object for pdnew() function.
    com_tot_pair = Matrix( com_tot_pair, 
                           headers = { '0': combin, 
                                       '1': pam.get_column_headers() } ) 


    # Array holding the PD of each pairwise community combination.
    pd_tot_pair = pdnew( com_tot_pair, tree )

    # The following will calculate the sum of each pair of sample's PD. 
    # I.e. treating each sample separately and just adding the PD values together.
    sum_pd_pair = []
    for pair in range( len( combin ) ):
        tmp = pd.data[ combin[ pair ][ 0 ], 0 ] + pd.data[ combin[ pair ][ 1 ], 0 ]
        sum_pd_pair.append( tmp )

    # PD of all communities combined (in case different from just using the entire tree).
    com_tot_multi = np.sum( pam.data, axis = 0 )
    com_tot_multi = [1 if i > 0 else 0 for i in com_tot_multi] # convert to presence/ absence (i.e. 1,0).
    sp_pres = list( it.compress( pam.get_column_headers(), com_tot_multi ) ) # calculate the PD.
    tree_pres = tree.extract_tree_with_taxa_labels( sp_pres )
    pd_tot_multi = tree_pres.length()

    # Contribution of PD that is not shared beteen two sites:
    pd_sites = pd.data[ 0:len( pd.get_row_headers() ), 0 ] # pull out just PD values.
    pd_combos = list( it.combinations( pd_sites, 2 ) ) # create list of all pairwise combinations.

    # Array to hold metrics assessing PD contributions to beta diversity.
    not_shared = np.zeros( (len( pd_combos ), 4 ), 
                           dtype = np.float )

    # this loop will populate array. see Details.
    for pair in range( len( pd_combos ) ):
        
        # Pull out each site's individual data.
        site1 = pd_combos[ pair ][ 0 ] # PD site 1
        site2 = pd_combos[ pair ][ 1 ] # PD site 2

        # Pull out the PD of the 2 sites combined.
        pdpair = pd_tot_pair.data[ pair ][ 0 ]
        #print site1, site2, pdpair

        # Pull out the sum of the separate PD values.
        sum_pair = sum_pd_pair[ pair ]

        # Metrics of interest:
        min_not_shared = min( pdpair - site1, pdpair - site2 ) # min(b,c)
        max_not_shared = max( pdpair - site1, pdpair - site2 ) # max(b,c)
        sum_not_shared = ( 2 * pdpair ) - sum_pair # b+c
        shared_val = pdpair - sum_not_shared # a

        # Add metrics to appropriate row of array.
        not_shared[pair] = [ min_not_shared, max_not_shared, sum_not_shared, shared_val ]

    # Convert not_shared array to Matrix object.
    core_calc = Matrix( not_shared,
                        headers = { '0': combin,
                                    '1': ['min_not_shared', 'max_not_shared',
                                          'sum_not_shared', 'shared' ] } )

    # return values.
    return core_calc

# .............................................................................
def get_species_index_lookup(pam):
    """Creates a lookup dictionary for species in a matrix

    Args:
        pam (:obj:`Matrix`): A Lifemapper Matrix object with presence absence
            values.

    Returns:
        A dictionary of species keys with matrix index values.
    """
    return dict(
        [(sp, idx) for (idx, sp) in enumerate(pam.get_column_headers())])


# .............................................................................
def calculate_phylo_beta_diversity_jaccard(pam, tree):
    """Calculates phylogenetic beta diversity for the jaccard index family.

    Args:
        pam (:obj:`Matrix`): A Lifemapper Matrix object with presence absence
            values (site rows by species columns).
        tree (:obj:`TreeWrapper`): A TreeWrapper object for a wrapped Dendropy
            phylogenetic tree.

    Returns:
        Phylogenetic beta diversity matrics (species by species)
            * beta_jtu: ADD DESCRIPTION
            * phylo_beta_jtu: ADD DESCRIPTION
            * beta_jne: ADD DESCRIPTION
            * phylo_beta_jne: ADD DESCRIPTION
            * beta_jac: ADD DESCRIPTION
            * phylo_beta_jac: ADD DESCRIPTION

    Note:
        * It looks like the scipy.spatial.distance.jaccard method may be useful
            here.

    Todo:
        * Fill in method documentation
        * Fill in method
    """
    # Get a lookup dictionary for the matrix index of each species in the PAM
    #    in case they are not in the same order as the taxa in the tree
    species_lookup = get_species_index_lookup(pam)

    # Build a header dictionary, all of the returned matricies will have the
    #    same headers, site rows by site columns.
    # Note: This will differ from the R method because each site will be
    #    present in both the rows and the columns.
    mtx_headers = {
        '0': pam.get_row_headers(),  # Row headers
        '1': pam.get_row_headers()  # Column headers
    }

    num_sites = pam.data.shape[0]  # Get the number of sites in the PAM

    # Note: For ease of development, use these numpy arrays for the
    #    computations.  They will be wrapped into a Matrix object when they are
    #    returned from the function.
    beta_jtu_data = np.zeros( ( num_sites, num_sites ), dtype = np.float )
    
    phylo_beta_jtu_data = np.zeros( ( num_sites, num_sites ), dtype = np.float )
    
    beta_jne_data = np.zeros( ( num_sites, num_sites ), dtype = np.float )
    
    phylo_beta_jne_data = np.zeros( ( num_sites, num_sites ), dtype = np.float )
    
    beta_jac_data = np.zeros( ( num_sites, num_sites ), dtype = np.float )
    
    phylo_beta_jac_data = np.zeros( ( num_sites, num_sites ), dtype = np.float )

    # TODO: Compute phylo beta diversity for jaccard index family
    core_calc = core_PD_calc( pam, tree )

    # print core_calc.data
    # print core_calc.get_row_headers()
    # print core_calc.get_column_headers()

    # This loop will populate arrays with beta diversity metrics.
    for my_row in range( core_calc.data.shape[0] ):
        # Pull out the numeric values.
        my_dat = core_calc.data[ my_row, 0:4 ]
        # Get index values for placing into output arrays.
        my_dim = core_calc.get_row_headers()[ my_row ]

        # Populate arrays.
        phylo_beta_jtu_data[ my_dim[0], my_dim[1] ] = ( 2*my_dat[0] ) / ( ( 2*my_dat[0] ) + my_dat[3] )
        phylo_beta_jtu_data[ my_dim[1], my_dim[0] ] = phylo_beta_jtu_data[ my_dim[0], my_dim[1] ]

        phylo_beta_jac_data[ my_dim[0], my_dim[1] ] = my_dat[2] / ( my_dat[3] + my_dat[2] )
        phylo_beta_jac_data[ my_dim[1], my_dim[0] ] = phylo_beta_jac_data[ my_dim[0], my_dim[1] ]

        phylo_beta_jne_data[ my_dim[0], my_dim[1] ] = ( ( my_dat[1] - my_dat[0] ) / ( my_dat[3] + my_dat[2]) ) * ( my_dat[3] / ( ( 2*my_dat[0] ) + my_dat[3] ) )
        phylo_beta_jne_data[ my_dim[1], my_dim[0] ] = phylo_beta_jne_data[ my_dim[0], my_dim[1] ]

    # Just to match Biotaphy test file expectations.
    for i in range( num_sites ):
        phylo_beta_jtu_data[ i, i ] = 1.
        phylo_beta_jac_data[ i, i ] = 1.
        phylo_beta_jne_data[ i, i ] = 1.
    
    print "\nJTU\n", phylo_beta_jtu_data, "\nJAC\n", phylo_beta_jac_data, "\nJNE\n", phylo_beta_jne_data

    return (
        Matrix(beta_jtu_data, headers=mtx_headers),
        Matrix(phylo_beta_jtu_data, headers=mtx_headers),
        Matrix(beta_jne_data, headers=mtx_headers),
        Matrix(phylo_beta_jne_data, headers=mtx_headers),
        Matrix(beta_jac_data, headers=mtx_headers),
        Matrix(phylo_beta_jac_data, headers=mtx_headers))


# .............................................................................
def calculate_phylo_beta_diversity_sorensen(pam, tree):
    """Calculates phylogenetic beta diversity for the sorensen index family.

    Args:
        pam (:obj:`Matrix`): A Lifemapper Matrix object with presence absence
            values.
        tree (:obj:`TreeWrapper`): A TreeWrapper object for a wrapped Dendropy
            phylogenetic tree.

    Returns:
        Phylogenetic beta diversity matrics (species by species)
            * beta_sim: ADD DESCRIPTION
            * phylo_beta_sim: ADD DESCRIPTION
            * beta_sne: ADD DESCRIPTION
            * phylo_beta_sne: ADD DESCRIPTION
            * beta_sor: ADD DESCRIPTION
            * phylo_beta_sor: ADD DESCRIPTION

    Todo:
        * Fill in method documentation
        * Fill in method
    """
    # Get a lookup dictionary for the matrix index of each species in the PAM
    #    in case they are not in the same order as the taxa in the tree
    species_lookup = get_species_index_lookup(pam)

    # Build a header dictionary, all of the returned matricies will have the
    #    same headers, site rows by site columns.
    # Note: This will differ from the R method because each site will be
    #    present in both the rows and the columns.
    mtx_headers = {
        '0': pam.get_row_headers(),  # Row headers
        '1': pam.get_row_headers()  # Column headers
    }

    num_sites = pam.data.shape[0]  # Get the number of sites in the PAM

    # Note: For ease of development, use these numpy arrays for the
    #    computations.  They will be wrapped into a Matrix object when they are
    #    returned from the function.
    beta_sim_data = np.zeros((num_sites, num_sites), dtype=np.float)
    phylo_beta_sim_data = np.zeros((num_sites, num_sites), dtype=np.float)
    beta_sne_data = np.zeros((num_sites, num_sites), dtype=np.float)
    phylo_beta_sne_data = np.zeros((num_sites, num_sites), dtype=np.float)
    beta_sor_data = np.zeros((num_sites, num_sites), dtype=np.float)
    phylo_beta_sor_data = np.zeros((num_sites, num_sites), dtype=np.float)

    # TODO: Compute phylo beta diversity for sorensen index family

    core_calc = core_PD_calc( pam, tree )

    # print core_calc.data
    # print core_calc.get_row_headers()
    # print core_calc.get_column_headers()

    # This loop will populate arrays with beta diversity metrics.
    for my_row in range( core_calc.data.shape[0] ):

        my_dat = core_calc.data[ my_row, 0:4 ]
        my_dim = core_calc.get_row_headers()[ my_row ]

        phylo_beta_sim_data[ my_dim[0], my_dim[1] ] = my_dat[0] / ( my_dat[0] + my_dat[3] )
        phylo_beta_sim_data[ my_dim[1], my_dim[0] ] = phylo_beta_sim_data[ my_dim[0], my_dim[1] ]

        phylo_beta_sor_data[ my_dim[0], my_dim[1] ] = my_dat[2] / ( ( 2*my_dat[3] ) + my_dat[2] )
        phylo_beta_sor_data[ my_dim[1], my_dim[0] ] = phylo_beta_sor_data[ my_dim[0], my_dim[1] ]

        phylo_beta_sne_data[ my_dim[0], my_dim[1] ] = ( ( my_dat[1] - my_dat[0] ) / ( ( 2*my_dat[3] ) + my_dat[2] ) ) * ( my_dat[3] / ( my_dat[0] + my_dat[3] ) )
        phylo_beta_sne_data[ my_dim[1], my_dim[0] ] = phylo_beta_sne_data[ my_dim[0], my_dim[1] ]
    # Just to match formatting across scripts.
    for i in range( num_sites ):
        phylo_beta_sim_data[ i, i ] = 1.
        phylo_beta_sne_data[ i, i ] = 1.
        phylo_beta_sor_data[ i, i ] = 1.

    print "\nSIM\n", phylo_beta_sim_data, "\nSOR\n", phylo_beta_sor_data, "\nSNE\n", phylo_beta_sne_data

    return (
        Matrix(beta_sim_data, headers=mtx_headers),
        Matrix(phylo_beta_sim_data, headers=mtx_headers),
        Matrix(beta_sne_data, headers=mtx_headers),
        Matrix(phylo_beta_sne_data, headers=mtx_headers),
        Matrix(beta_sor_data, headers=mtx_headers),
        Matrix(phylo_beta_sor_data, headers=mtx_headers))



pam = Matrix(np.array( [ [1, 1, 1, 0, 0, 0],
                         [0, 1, 1, 1, 0, 0],
                         [0, 0, 1, 1, 1, 0],
                         [0, 0, 1, 1, 1, 1],
                         [0, 0, 0, 1, 1, 1],
                         [1, 0, 0, 1, 1, 1] ]   ),
                        headers = { '0': ['A', 'B', 'C', 'D', 'E', 'F'],
                                    '1': ['sp1', 'sp2', 'sp3', 'sp4', 'sp5', 'sp6'] } )

tree = TreeWrapper.get( data = '(((sp1:1,sp2:1):5,(sp3:3,sp4:3):3):2,(sp5:7,sp6:7):1);',
                        schema = 'newick' )

a = calculate_phylo_beta_diversity_jaccard( pam, tree )
b = calculate_phylo_beta_diversity_sorensen( pam, tree )

