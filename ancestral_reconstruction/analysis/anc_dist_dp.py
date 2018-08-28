"""
@summary: Module containing code for calculating ancestral distributions
"""
import argparse
from datetime import datetime
import dendropy as dp
import math
from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import stats
import scipy.linalg as la
from scipy.stats import gaussian_kde as kde
import sys

import aln_reader

rcParams.update({'figure.autolayout': True})
plt.style.use('fivethirtyeight')


# .............................................................................
def generate_argparser():
    parser = argparse.ArgumentParser(
                prog="anc_distr_rec",
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-t", "--treefile", type=str, nargs=1, required=True,
                        help=("Input tree in Newick format."))
    parser.add_argument("-d", "--datafile", type=open, nargs=1, required=True,
                        help=("Datafile in modified Phylip with columns space separated.\
                                If not present, then data will be simulated."))
    parser.add_argument("--precut", action="store_true",
                        help=("Are the data already in predetermined \
                                categories (otherwise, it is expected that the\
                                 data are independent points)"))
    parser.add_argument("-o", "--outdir", type=os.path.abspath, nargs=1,
                        required=True, help=("Output directory."))
    parser.add_argument("-c", "--ncats", type=int, default=30,
                        help=("The number of categories (how much to split \
                                it up)."))
    parser.add_argument("--sumtoone", action="store_true", default=True,
                        help=("Should we make the distributions sum to one?"))
    parser.add_argument("--printtree", action="store_true", required=False,
                        help=("Show the tree (required ete3)."))
    parser.add_argument("--printrates", action="store_true",
                        help=("Print the rates in the outdir (requires \
                                matplotlib)?"))
    parser.add_argument("--printplots", action="store_true", default=True,
                        help=("Print the rates in the outdir (requires \
                                matplotlib)?"))
    return parser


# .............................................................................
def calc_cont_anc_states_bin(tree, char):
    """
    @summary: Calculate continuous ancestral states
    """
    df = 0
    nodenum = {}
    count = 0
    for k in tree.postorder_edge_iter():
        i = k.head_node
        if len(i.child_nodes()) == 0:
            i.data['val'] = float(i.data['cont_values'][char])
            i.data['valse'] = float(i.data['cont_values'][char])
        else:
            nodenum[i] = count
            count += 1
            df += 1
            i.data['val'] = 0.
            i.data['valse'] = 0
    df -= 1
    # compute the mlest of the root
    fullMcp = np.zeros((df+1, df+1))
    fullVcp = np.zeros(df+1)
    count = 0
    for k in tree.postorder_edge_iter():
        i = k.head_node
        if len(i.child_nodes()) > 0:
            nni = nodenum[i]
            for j in i.child_nodes():
                tbl = 2./j.edge_length
                fullMcp[nni][nni] += tbl
                if len(j.child_nodes()) == 0:
                    fullVcp[nni] += (j.data['val'] * tbl)
                else:
                    nnj = nodenum[j]
                    fullMcp[nni][nnj] -= tbl
                    fullMcp[nnj][nni] -= tbl
                    fullMcp[nnj][nnj] += tbl
            count += 1
    # NOTE: these two cholesky steps are the slowest bits
    b = la.cho_factor(fullMcp)
    # these are the ML estimates for the ancestral states
    mle = la.cho_solve(b, fullVcp)
    sos = 0
    for k in tree.postorder_edge_iter():
        i = k.head_node
        if len(i.child_nodes()) > 0:
            i.data['val'] = mle[nodenum[i]]
            i.data['cont_values'].append(mle[nodenum[i]])
            # print(i.data['val'])
            # i.label = str(mle[nodenum[i]])
            # print(i.get_newick_repr(False), mle[nodenum[i]])
            for j in i.child_nodes():
                temp = (i.data['val'] - j.data['val'])
                sos += temp*temp / j.edge_length
    # calcSE
    # need to change this to the rohlf 2001 standard errors
    for k in tree.postorder_edge_iter():
        i = k.head_node
        if len(i.child_nodes()) > 0:
            qpq = fullMcp[nodenum[i]][nodenum[i]]
            tm1 = np.delete(fullMcp, (nodenum[i]), axis=0)
            tm = np.delete(tm1, (nodenum[i]), axis=1)
            b = la.cho_factor(tm)
            sol = la.cho_solve(b, tm1[:, nodenum[i]])
            tempse = qpq - np.inner(tm1[:, nodenum[i]], sol)
            i.data['valse'] = math.sqrt(2*sos/(df*tempse))
            i.data['cont_values_se'].append(i.data['valse'])
            # print(i.data['valse'])


# .............................................................................
def match_tips_and_cont_values(tree, seqs):
    ln = {}
    sq = {}
    for i in tree:
        i.data = {}
        if len(i.child_nodes()) == 0:
            ln[i.taxon.label] = i
    for i in seqs:
        sq[i.name] = i
    for i in ln:
        try:
            ln[i].data['cont_values'] = sq[i].cont_values
            ln[i].data['orig_values'] = sq[i].orig_values
        except:
            print("Can't find {} in cont_values".format(i.label))
            return False


# .............................................................................
def scale_to_one(li):
    sc = float(1)/sum(li)
    nl = []
    for i in li:
        nl.append(i*sc)
    return nl, sc


# .............................................................................
def scale_and_se(tree, args):
    for k in tree.postorder_edge_iter():
        i = k.head_node
        if len(i.child_nodes()) > 0:
            if args.sumtoone:
                i.data['cont_values'], scale = scale_to_one(
                                                        i.data['cont_values'])
                i.data['cont_values_se'] = [j*scale for j in i.data[
                                                            'cont_values_se']]
            i.data['cont_values_low'] = [max(j-k, 0) for j, k in zip(
                            i.data['cont_values'], i.data['cont_values_se'])]
            i.data['cont_values_high'] = [max(j+k, 0) for j, k in zip(
                            i.data['cont_values'], i.data['cont_values_se'])]


# .............................................................................
def print_rates(rates, x_grid, outd):
    plt.plot(x_grid, rates)
    plt.ylabel('estimated rate')
    plt.xlabel('state space')
    plt.savefig(outd+"rate_estimate.png")
    # plt.show()


# .............................................................................
def calculate_densities(seqs, low, high, args):
    x_grid = np.linspace(low, high, args.ncats)
    for i in seqs:
        if not args.precut or args.datafile is None:
            density = kde(i.cont_values)
            density.covariance_factor = lambda: .25
            density._compute_covariance()
            kdepdf = density.evaluate(x_grid)
            i.orig_values = i.cont_values
            if args.sumtoone:
                i.cont_values, scale = scale_to_one(i.cont_values)
            i.cont_values = kdepdf
        else:
            i.orig_values = i.cont_values
            if args.sumtoone:
                i.cont_values, scale = scale_to_one(i.cont_values)
    return x_grid


# .............................................................................
def print_plots_to_file(tree, x_grid, args, outd):
    print("creating figures now")
    # Construct the png plot figures of the tips and internal nodes
    ndcount = 0
    for k in tree.postorder_edge_iter():
        i = k.head_node
        plt.figure(figsize=(6, 4))
        plt.plot(x_grid, i.data['cont_values'])
        # should output the i.data['cont_values'] for each node here
        if len(i.child_nodes()) > 0:
            i.label = "nd"+str(ndcount)
            ndcount += 1
        # else:
        #     if args.precut == False:
        #         plt.hist(i.data['orig_values'],normed=1,
        #         histtype='stepfilled', alpha=0.25)
        plt.fill_between(x_grid, 0, i.data['cont_values'], alpha=0.05)
        if len(i.child_nodes()) > 0:
            plt.plot(x_grid, i.data['cont_values_high'], '--', alpha=0.55,)
            plt.plot(x_grid, i.data['cont_values_low'], '--', alpha=0.55,)
        plt.grid(True)
        if i.label is None:
            plt.savefig(outd+str(i.taxon.label)+'.png')
        else:
            plt.savefig(outd+str(i.label)+'.png')
        plt.close()


# .............................................................................
def main():
    arguments = sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(arguments)

    start = datetime.now()
    tree = dp.Tree.get(path=args.treefile[0], schema="newick")
    outd = args.outdir[0]
    if outd[-1] != "/":
        outd += "/"
    print("out dir: {}".format(outd))
    if not os.path.isdir(outd):
        os.makedirs(outd)

    # read data file if present
    low = None
    high = None
    seqs = None
    cats = []
    if args.datafile:
        seqs = aln_reader.read_table_cont_file(args.datafile[0])
        if args.precut:
            ncuts = len(seqs[0].cont_values)
            args.ncats = ncuts
            cats = list(range(0, args.ncats))
            low = 0.
            high = args.ncats
        else:
            exit(0)
    else:
        exit(0)

    # this calculates the categories if not precut
    if not args.precut:
        for i in range(args.ncats+1):
            cats.append(low + (i*(float(high-low)/args.ncats)))
    print("low range: {}".format(low))
    print("high range:{}".format(high))
    print("cats:{}".format(cats))

    # Calculate the kernel densities
    print("calculating densisties")
    x_grid = calculate_densities(seqs, low, high, args)
    print("finished calculating densities")
    match_tips_and_cont_values(tree, seqs)

    # prepare the tree
    for k in tree.postorder_edge_iter():
        i = k.head_node
        if len(i.child_nodes()) != 0:
            i.data['cont_values'] = []
            i.data['cont_values_se'] = []

    # Conduct the analyses
    rates = []
    print("calculating anc states")
    for i in range(args.ncats):
        print(" cat:"+str(i))
        calc_cont_anc_states_bin(tree, i)
        # estrate = sigsqML(tree)
        # rates.append(estrate)
    print("\nfinished calculation")

    # plot the estimated rates
    if args.printrates:
        print_rates(rates, x_grid, outd)

    # standard error and scale to one
    scale_and_se(tree, args)

    end = datetime.now()
    print("runtime for analysis (H:M:S): "+str(end-start))

    if args.printplots:
        print_plots_to_file(tree, x_grid, args, outd)
        ts = tree.as_string("newick")
        outf = open(outd+"treefile_labeled.tre", "w")
        outf.write(ts)
        outf.close()

    # Plotting the tree and distribution using ETE
    if args.printtree and args.printplots:
        print_tree_to_file(tree, outd)


# .............................................................................
if __name__ == "__main__":
    main()
