import sys
import pre_process

if __name__ == '__main__':
    #debug: args = str(sys.argv)
    d = int(sys.argv[1])
    l = int(sys.argv[2])
    min_sup = int(sys.argv[3])
    query = sys.argv[4:]

    #debug: print (query)
    dataSet = {}  # {window: genome id}
    dataSet = pre_process.start(d, query)

    print(dataSet)
    # 2.build tree ( # build tree from windows)
    # 3.return ans(get ans from tree) (# return hitchhikers ??)
