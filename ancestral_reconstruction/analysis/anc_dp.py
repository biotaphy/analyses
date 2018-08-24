"""
@summary: Module containing code for calculating ancestral states
"""
import dendropy as dp
import math
import numpy as np
import scipy.linalg as la
import sys

"""
sq_change which produces the same results as ML without SE
"""
def calc_cont_anc_states(tree):
    """
    @summary: Calculate continuous ancestral states for tree nodes
    """
    df = 0
    nodenum = {}
    count = 0
    for k in tree.postorder_edge_iter():
        i = k.head_node
        if len(i.child_nodes()) == 0:
            i.data['val'] = float(i.data['cont_values'][0])
            i.data['valse'] = float(i.data['cont_values'][0])
        else:
            nodenum[i] = count
            count += 1
            df += 1
            i.data['val'] = 0.
            i.data['valse'] = 0.
    df -= 1
    #compute the mlest of the root
    fullMcp = np.zeros((df+1,df+1))
    fullVcp = np.zeros(df+1)
    count = 0
    for k in tree.postorder_edge_iter():
        i = k.head_node
        if len(i.child_nodes()) != 0:
            nni = nodenum[i]
            for j in i.child_nodes():
                tbl = 2./j.edge_length
                fullMcp[nni][nni] += tbl;
                if len(j.child_nodes()) == 0:
                    fullVcp[nni] += (j.data['val'] * tbl)
                else:
                    nnj = nodenum[j]
                    fullMcp[nni][nnj] -= tbl;
                    fullMcp[nnj][nni] -= tbl;
                    fullMcp[nnj][nnj] += tbl;
            count += 1
    b = la.cho_factor(fullMcp)
    #these are the ML estimates for the ancestral states
    mle = la.cho_solve(b,fullVcp)
    sos = 0
    for k in tree.postorder_edge_iter():
        i = k.head_node
        if len(i.child_nodes()) != 0:
            i.data['val'] = mle[nodenum[i]]
            #print i.data['val']
            i.label = str(mle[nodenum[i]])
            for j in i.child_nodes():
                temp = (i.data['val'] - j.data['val'])
                sos += temp*temp / j.edge_length
    #print "Square Length: ",sos
    #calcSE
    """
    for i in tree.iternodes(order="postorder"):
        if i.istip == False:
            qpq = fullMcp[nodenum[i]][nodenum[i]]
            tm1 = np.delete(fullMcp,(nodenum[i]),axis=0)
            tm = np.delete(tm1,(nodenum[i]),axis=1)
            b = cho_factor(tm)
            sol = cho_solve(b,tm1[:,nodenum[i]])
            tempse = qpq - np.inner(tm1[:,nodenum[i]],sol)
            i.data['valse'] = math.sqrt(2*sos/(df*tempse))
    """
    return 0



def match_tips_and_cont_values(tree,seqs):
    for i in tree:
        i.data = {}
        if len(i.child_nodes()) == 0:
            test = False
            for j in seqs:
                if i.taxon.label == j.name:
                    test = True
                    i.data['cont_values'] = j.cont_values
                    break
            if test == False:
                print "can't find "+i.taxon.label+" in cont_values"
                return False

# .............................................................................
def calculate_continuous_ancestral_states(tree, sequences):
    """
    @summary: Calculates continuous ancestral states for a tree and a list of
                 sequences
    @note: sq_change which produces the same results as ML without SE
    """
    match_tips_and_cont_values(tree, sequences)
    df = 0
    nodenum = {}
    count = 0
    for k in tree.postorder_edge_iter():
        i = k.head_node
        if len(i.child_nodes()) == 0:
            i.data['val'] = float(i.data['cont_values'][0])
            i.data['valse'] = float(i.data['cont_values'][0])
        else:
            nodenum[i] = count
            count += 1
            df += 1
            i.data['val'] = 0.
            i.data['valse'] = 0.
    df -= 1
    #compute the mlest of the root
    fullMcp = np.zeros((df+1,df+1))
    fullVcp = np.zeros(df+1)
    count = 0
    for k in tree.postorder_edge_iter():
        i = k.head_node
        if len(i.child_nodes()) != 0:
            nni = nodenum[i]
            for j in i.child_nodes():
                tbl = 2./j.edge_length
                fullMcp[nni][nni] += tbl;
                if len(j.child_nodes()) == 0:
                    fullVcp[nni] += (j.data['val'] * tbl)
                else:
                    nnj = nodenum[j]
                    fullMcp[nni][nnj] -= tbl;
                    fullMcp[nnj][nni] -= tbl;
                    fullMcp[nnj][nnj] += tbl;
            count += 1
    b = la.cho_factor(fullMcp)
    #these are the ML estimates for the ancestral states
    mle = la.cho_solve(b,fullVcp)
    sos = 0
    for k in tree.postorder_edge_iter():
        i = k.head_node
        if len(i.child_nodes()) != 0:
            i.data['val'] = mle[nodenum[i]]
            #print i.data['val']
            i.label = str(mle[nodenum[i]])
            for j in i.child_nodes():
                temp = (i.data['val'] - j.data['val'])
                sos += temp*temp / j.edge_length
    #print "Square Length: ",sos
    #calcSE
    """
    for i in tree.iternodes(order="postorder"):
        if i.istip == False:
            qpq = fullMcp[nodenum[i]][nodenum[i]]
            tm1 = np.delete(fullMcp,(nodenum[i]),axis=0)
            tm = np.delete(tm1,(nodenum[i]),axis=1)
            b = cho_factor(tm)
            sol = cho_solve(b,tm1[:,nodenum[i]])
            tempse = qpq - np.inner(tm1[:,nodenum[i]],sol)
            i.data['valse'] = math.sqrt(2*sos/(df*tempse))
    """
    return 0



#if __name__ == "__main__":
#    if len(sys.argv) != 3:
#        print "python "+sys.argv[0]+" newick.tre dataasphy"
#        sys.exit(0)
#    tree = dp.Tree.get(path=sys.argv[1],schema="newick")    
#    seqs = aln_reader.read_phylip_cont_file(open(sys.argv[2],"r"))
#    match_tips_and_cont_values(tree,seqs)
#    sqch = calc_cont_anc_states(tree)
#    outfile = open("contanc_dp.tre","w")
#    outfile.write(tree.as_string(schema="newick"))
#    outfile.close()
