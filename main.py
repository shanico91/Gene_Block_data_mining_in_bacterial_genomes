import sys
import pre_process
import tree_node
import tree_funcs

if __name__ == '__main__':
    #debug: args = str(sys.argv)
    d = int(sys.argv[1])
    l = int(sys.argv[2])
    min_sup = int(sys.argv[3])
    query = sys.argv[4:]

    # debug: print (query)
    data_set = {}  # {window: genome id}
    data_set = pre_process.start(d, query)
    print(data_set)

    #tree, header_table = tree_funcs.create_tree(data_set, min_sup)
    #tree.disp()
    # 2.build tree ( # build tree from windows)
    # 3.return ans(get ans from tree) (# return hitchhikers ??)
