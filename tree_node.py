# variables:
# name of the node, a count
# node_link used to link similar items
# parent variable used to refer to the parent of the node in the tree
# node contains an empty dictionary for the children in the node


class TreeNode:

    def __init__(self, name_value, num_occur, parent_node):
        self.name = name_value  # gene cog number
        self.count = num_occur  # in how many genomes this gene appeared in a valid window
        self.node_link = None  # points to the next node with the same cog number
        self.parent = parent_node
        self.children = {}
        self.last_update = 0  # the last genome that updated this node in the creation of the tree
        if parent_node is None:  # will be used to apply the l (min size of report) constriction
            self.depth = 0
        else:
            self.depth = parent_node.depth + 1

    # increments the count variable with a given amount

    def inc(self, num_occur):
        self.count += num_occur

    # display tree in text. Useful for debugging

    def disp(self, ind=1):
        print('  ' * ind, self.name, ' ', self.count)
        for child in self.children.values():
            child.disp(ind + 1)
