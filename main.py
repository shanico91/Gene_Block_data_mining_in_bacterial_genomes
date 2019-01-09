import sys
import pre_process
import tree_node
import tree_funcs

if __name__ == '__main__':
    # debug: args = str(sys.argv)
    d = int(sys.argv[1])
    length = int(sys.argv[2])  # l
    min_sup = int(sys.argv[3])
    query = sys.argv[4:]

    # debug: print (query)
    data_set = {}  # {window: genome id}
    data_set = pre_process.start(d, query, length)
    # print(data_set)

    tree = tree_funcs.create_tree(data_set, length-len(query), min_sup)
    freq_paths = []
    tree_funcs.dfs(tree, freq_paths, [])
    print(freq_paths)
    # tree.disp()