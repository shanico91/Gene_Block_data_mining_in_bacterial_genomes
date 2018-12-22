import tree_node


# createTree(), takes the dataset and the minimum support as arguments and builds the FP-tree.
# This makes two passes through the dataset.
# The first pass goes through everything in the dataset and counts the frequency of each item.
# These are stored in the header table.
def create_tree(data_set, min_sup=1):  # create FP-tree from dataset but don't mine
    header_table = {}  # {gene: no. of genomes it appears in}
    # go over dataSet twice:

    # first pass counts frequency of occurrence
    for genome in data_set:  # genome id
        local_set = set()
        for window in data_set[genome]:
            for gene in window:
                local_set.add(gene)
        for gene in local_set:
            header_table[gene] = header_table.get(gene, 0) + 1

    for gene in list(header_table):  # remove items not meeting minSup
        if header_table[gene] < min_sup:
            del (header_table[gene])
    freq_gene_set = set(header_table.keys())
    # print 'freq_gene_set: ',freq_gene_set

    if len(freq_gene_set) == 0:
        return None, None  # if no genes meet min support -->get out

    for gene in header_table:
        header_table[gene] = [header_table[gene], None]  # reformat header_table to use Node link
    # print 'header_table: ',header_table

    ret_tree = tree_node.TreeNode('Null Set', 1, None)  # create tree

    for genome, windows in data_set.items():
        for window in windows:
            local_d = {}  # {gene: no of genomes it's in}
            for gene in window:
                if gene in freq_gene_set:
                    local_d[gene] = header_table[gene][0]
                    # TODO: this makes sure that each gene appears only once in ordered items- correct?
            if len(local_d) > 0:  # there are frequent genes in this window
                ordered_items = [v[0] for v in sorted(local_d.items(), key=lambda p: p[1], reverse=True)]
                # (ordered_items is the window that is ordered in a descending order based on
                # frequency of each gene in the genomes)
                update_tree(ordered_items, ret_tree, header_table, genome)  # populate tree with ordered freq itemset

    return ret_tree, header_table  # return tree and header table


# updateTree() grow the Fp-tree with a window of genes.
def update_tree(window, in_tree, header_table, genome_num):
    if window[0] in in_tree.children:  # check if orderedItems[0] in retTree.children
        if in_tree.children[window[0]].last_update != genome_num:
            in_tree.children[window[0]].inc(1)  # increase count by one
            in_tree.children[window[0]].last_update = genome_num
    else:  # add items[0] to in_tree.children
        in_tree.children[window[0]] = tree_node.TreeNode(window[0], 1, in_tree)
        in_tree.children[window[0]].last_update = genome_num
        # update header table
        if header_table[window[0]][1] is None:
            header_table[window[0]][1] = in_tree.children[window[0]]
        else:
            update_header(header_table[window[0]][1], in_tree.children[window[0]])
    if len(window) > 1:  # call updateTree() with remaining ordered items
        update_tree(window[1::], in_tree.children[window[0]], header_table, genome_num)


# updateHeader() makes sure the node links points to every instance of the this item on the tree
def update_header(node_to_test, target_node):  # this version does not use recursion
    while node_to_test.node_link is not None:  # Do not use recursion to traverse a linked list!
        node_to_test = node_to_test.node_link
    node_to_test.node_link = target_node


# ascendTree(), which ascends the tree, collecting the names of items it encounters
def ascend_tree(leaf_node, prefix_path):  # ascends from leaf node to root
    if leaf_node.parent is not None:
        prefix_path.append(leaf_node.name)
        ascend_tree(leaf_node.parent, prefix_path)


# The findPrefixPath() function iterates through the linked list until it hits the end.
# For each item it encounters, it calls ascendTree().
# This list is returned and added to the conditional pattern base dictionary called condPats.
def find_prefix_path(base_pat, node):  # treeNode comes from header table
    # TODO: figure out what is the use of base_pat
    cond_pats = {}
    while node is not None:
        prefix_path = []
        ascend_tree(node, prefix_path)
        if len(prefix_path) > 1:
            cond_pats[frozenset(prefix_path[1:])] = node.count
        node = node.nodeLink
    return cond_pats
