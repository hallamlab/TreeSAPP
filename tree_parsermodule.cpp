#include <Python.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <stack>

using namespace std;

static PyObject *read_the_reference_tree(PyObject *self, PyObject *args);
static PyObject *get_parents_and_children(PyObject *self, PyObject *args);
static PyObject *build_subtrees_newick(PyObject *self, PyObject *args);
static PyObject *lowest_common_ancestor(PyObject *self, PyObject *args);
char *get_node_relationships(char *tree_string);
char *split_tree_string(char *tree_string);

static char module_docstring[] =
        "This module provides an interface for parsing Newick formatted trees using C from within MLTreeMap";
static char read_the_reference_tree_docstring[] =
        "Reads the labelled_tree_file and reformats it for downstream interpretation";
static char get_parents_and_children_docstring[] = 
        "Stores the input tree as a binary search tree before recursively finding the children and parent of each node";
static char build_subtrees_newick_docstring[] =
        "Reads the labelled, rooted tree and returns all subtrees in the tree";
static char lowest_common_ancestor_docstring[] =
        "Calculate lowest common ancestor for a set of nodes in a tree";

static PyMethodDef module_methods[] = {
        {"_read_the_reference_tree",
        read_the_reference_tree,
        METH_VARARGS,
        read_the_reference_tree_docstring},
        {"_get_parents_and_children",
        get_parents_and_children,
        METH_VARARGS,
        get_parents_and_children_docstring},
        {"_build_subtrees_newick",
        build_subtrees_newick,
        METH_VARARGS,
        build_subtrees_newick_docstring},
        {"_lowest_common_ancestor",
        lowest_common_ancestor,
        METH_VARARGS,
        lowest_common_ancestor_docstring},
        {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_tree_parser(void) {
    PyObject *m = Py_InitModule3("_tree_parser", module_methods, module_docstring);
    if (m == NULL)
        return;
}

/**
  A tree node
  */
struct Link {
    long key;
    long left;
    long right;
    Link* next;
    Link* previous;
};

struct TreeNode {
    long key;
    TreeNode* left;
    TreeNode* right;
};

struct CharLink {
    char* subtree;
    CharLink* next;
};


/**
 * Inserts a new Link (with key=newKey) at the head of the linked_list.
 */
void prepend_link(Link*& head, long newKey) {
    Link * curr = new Link;
    curr->key  = newKey;
    curr->next = head;
    head = curr;
}


TreeNode *& find(int key, TreeNode *& r) {
    if (r == NULL) return r;
    if (key < r->key)
        return find(key, r->left);
    if (key > r->key)
        return find(key, r->right);
    return r;
}


TreeNode* create_node(long key, TreeNode* l = NULL, TreeNode* r = NULL ) {
    TreeNode* curr = new TreeNode;
    curr->key = key;
    curr->left = l;
    curr->right = r;
    return curr;
}


void add_subtree(CharLink*& head, char*& new_subtree) {
    CharLink * curr = new CharLink;
    curr->subtree = new_subtree;
    curr->next = head;
    head = curr;
}


void deleteSubtreeList(CharLink*& head) {
    if ( head != NULL ) {
        deleteSubtreeList( head->next );
        free(head->subtree);
        delete head;
        head = NULL;
    }
}


void deleteList(Link*& head) {
    if ( head != NULL ) {
        deleteList( head->next );
        delete head;
        head = NULL;
    }
}


/**
 * Deletes all nodes in the tree rooted at root and sets root to NULL.
 */
void deleteTree( TreeNode*& root ) {
    if ( root != NULL ) {
        deleteTree( root->left );
        deleteTree( root->right );
        delete root;
        root = NULL;
    }
}


void print_list(Link* head) {
    std::cout << std::endl;
    for (Link* curr = head; curr != NULL; curr = curr->next){
        printf("key: %ld\tleft: %ld\tright: %ld",
            curr->key, curr->left, curr->right);
        if (curr->next != NULL) std::cout << "\n";
    }
    std::cout << std::endl;
}


/**
 * Prints out the tree sideways.
 */
void print_tree( TreeNode* root, int d = 0 ) {
    if ( root == NULL ) return;
    print_tree( root->right, d+1 );
    std::cout << std::setw( 3 * d ) << ""; // output 3 * d spaces
    std::cout << root->key << std::endl;
    print_tree( root->left, d+1 );
}

char const *float_chars = "0.123456789";
char const *real_number_chars = "-0123456789";

/*
 Returns 1 if character sub is a substring of super, 0 otherwise
 */
int is_char_substr(char sub, const char * super) {
    int x = 0;
    while(super[x]) {
        if (super[x] == sub)
            return 1;
        x++;
    }
    return 0;
}

/*
 Returns the length of the char_array as int
 */
int get_char_array_length(char * char_array) {
    int x = 0;
    while (char_array[x])
        x++;
    return x;
}

/*
 * Append a character array (source) onto another character array (dest)
 * start is the position to continue appending on dest
 */
int append_char_array(int start, char * source, char *&dest) {
    int i = 0;
    while (source[i])
        dest[start++] = source[i++];
    return start;
}

char reverse_char_array(char * char_array, char *&flipped, int first, int last) {
    if (last == -1)
        return '\0';
    flipped[first] = reverse_char_array(char_array, flipped, first+1, last-1);
    return char_array[last];
}

/*
 * param comma_separated_string: a character array with commas
 * return: a char** the elements stored as individual char arrays
 */
char **csv_to_list(char * comma_separated_string) {
    int i = 0;
    int j = 0;
    int k = 0;
    char** char_list = (char**) malloc (10 * sizeof(char*));

    char_list[j] = (char*) malloc(10);
    while (comma_separated_string[i]) {
        if (comma_separated_string[i] == ',') {
            char_list[j][k] = '\0';
            i++; j++; k = 0;
            char_list[j] = (char*) malloc(10);
        }
        char_list[j][k] = comma_separated_string[i];
        i++; k++;
    }
    char_list[j][k] = '\0';
    char_list[++j] = '\0';
    return char_list;
}


static PyObject *read_the_reference_tree(PyObject *self, PyObject *args) {
    char* reference_tree_file;
    if (!PyArg_ParseTuple(args, "s", &reference_tree_file)) {
        return NULL;
    }

    FILE *reference_tree = fopen(reference_tree_file, "r");
    if (reference_tree == 0) {
        printf ("The reference tree file %s could not be opened for reading!\n", reference_tree_file);
        exit(0);
    }

    int _MAX = 1000;
    int tree_len = 0;
    int count = 2;
    int x = 0;
    char *tree_string = (char *) malloc ( _MAX * sizeof(char));
    char count_char[10];
    char c = fgetc(reference_tree);
    while ( c != EOF ) {
        if (tree_len <= _MAX && tree_len >= _MAX - 100) {
            _MAX = _MAX + 1000;
            tree_string = (char *) realloc(tree_string, _MAX * sizeof(char));
        }
        if (c == ')') {
            tree_string[tree_len] = c;
            tree_len++;
            tree_string[tree_len] = '-';
            tree_len++;

            sprintf(count_char, "%d", count);
            count++;
            x = 0;
            while (count_char[x]) {
                c = count_char[x];
                tree_string[tree_len] = c;
                tree_len++;
                x++;
            }
            c = fgetc(reference_tree);
        }
        if (c == ':') {
            c = fgetc(reference_tree);
            // while c is a substring of float_chars, continue reading in characters
            while (is_char_substr(c, float_chars) == 1)
                c = fgetc(reference_tree);
        }
        else {
            tree_string[tree_len] = c;
            tree_len++;
            c = fgetc(reference_tree);
        }

    }
    fclose(reference_tree);
    // Now remove the last node
    while ( c != ')') {
        tree_string[tree_len + 1] = '\0';
        tree_string[tree_len] = ';';
        tree_len--;
        c = tree_string[tree_len];
    }

    return Py_BuildValue("s", tree_string);
};


void get_previous_node(char *&parsed_tree_string, int &end, char *&previous) {
    char reversed[10];
    char c = parsed_tree_string[end];
    int i = 0;
    // Skip through brackets and commas to the end of the previous node
    while (is_char_substr(c, real_number_chars) == 0) {
        parsed_tree_string[end--] = '\0';
        c = parsed_tree_string[end];
    }
    while (is_char_substr(c, real_number_chars) == 1) {
        reversed[i] = c;
        parsed_tree_string[end--] = '\0';
        c = parsed_tree_string[end];
        i++;
    }
    reversed[i] = '\0';
    reverse_char_array(reversed, previous, 0, i);
    return;
}


void load_linked_list(char * tree_string, Link *&head) {
    char c;
    int pos = 0;
    int i = 0;
    int newKey = -1;
    int retrace_pos = 0;
    int x;
    char curr[10];
    char* right = (char*) malloc(20);
    char* left = (char*) malloc(20);

    int tree_len = get_char_array_length(tree_string);
    char* parsed_tree_string = (char*) malloc(tree_len);
    for (x = 0; x < tree_len; x++)
        parsed_tree_string[x] = '\0';

    while (tree_string[pos]) {
        c = tree_string[pos];
        parsed_tree_string[retrace_pos++] = c;
        if (c == ')') {
            // load the next node as curr
            c = tree_string[pos+1];
            i = 0;
            // Overwrite curr
            for (x = 0; x < 10; x++)
                curr[x] = '\0';
            while (is_char_substr(c, real_number_chars) == 1) {
                curr[i++] = c;
                pos++;
                c = tree_string[pos+1];
            }
            curr[i] = '\0';
            newKey = atoi(curr);
            if (newKey == 0)
                newKey = -1;
            prepend_link(head, newKey);

            // load the previous 2 nodes as children and remove these from the string
            get_previous_node(parsed_tree_string, retrace_pos, right);
            head->right = atoi(right);
            for (x = 0; x < 10; x++)
                right[x] = '\0';
            get_previous_node(parsed_tree_string, retrace_pos, left);
            head->left = atoi(left);
            for (x = 0; x < 10; x++)
                left[x] = '\0';

            // add the current node to the parsed_tree_string
            i = 0;
            while (curr[i])
                parsed_tree_string[retrace_pos++] = curr[i++];
        }
        pos++;
    }
    free(right);
    free(left);
    free(parsed_tree_string);
}


TreeNode* load_tree_from_list(Link* head, TreeNode*& root, std::stack<TreeNode*>& merge) {
    TreeNode* previous = NULL;
    if (head == NULL) {
        return previous;
    }
    previous = load_tree_from_list(head->next, root, merge);
    root = create_node(head->key);

    if (previous) {
        // If a key is equal to the previous key, then the previous key is a child
        if (head->left == previous->key)
            root->left = previous;
        if (head->right == previous->key)
            root->right = previous;
        // If neither children are from the previous node, then create new nodes
        if (head->left != previous->key)
            root->left = create_node(head->left);
        if (head->right != previous->key)
            root->right = create_node(head->right);
        // If a child is equal to a node from a long time ago, in a galaxy far far away...
        if (!merge.empty() && head->left == merge.top()->key) {
            root->left = merge.top();
            merge.pop();
        }
        if (!merge.empty() && head->right == merge.top()->key) {
            root->right = merge.top();
            merge.pop();
        }
        // If neither the left or right are equal to the previous key, store it as merge
        if (head->right != previous->key && head->left != previous->key) {
            merge.push(previous);
        }
    }
    else {
        root->right = create_node(head->right);
        root->left = create_node(head->left);
    }
    return root;
}


int get_children_of_nodes(Link * head, char *&children) {
    char * buffer = (char*) malloc (20);
    int _MAX = 1000;
    children = (char *) malloc (_MAX * sizeof(char));
    int x = 0;
    for (Link* curr = head; curr != NULL; curr = curr->next){
        if (x <= _MAX && x >= _MAX - 100) {
            _MAX = _MAX + 1000;
            children = (char *) realloc (children, _MAX * sizeof(char));
        }
        sprintf(buffer, "%ld", curr->key);
        x = append_char_array(x, buffer, children);
        children[x++] = '=';
        sprintf(buffer, "%ld", curr->left);
        x = append_char_array(x, buffer, children);
        children[x++] = ',';
        sprintf(buffer, "%ld", curr->right);
        x = append_char_array(x, buffer, children);

        if (curr->next != NULL) children[x++] = ';';
    }
    children[x++] = '\n';
    children[x] = '\0';

    free(buffer);

    return x;
}


int get_parents_of_nodes(Link * head, char *&parents) {
    char parent_string[100];
    int _MAX = 1000;
    parents = (char *) malloc (_MAX * sizeof(char));
    int x;
    for (x = 0; x < _MAX; x++)
        parents[x] = '\0';
    for (x = 0; x < 100; x++)
        parent_string[x] = '\0';
    x = 0;
    for (Link* curr = head; curr != NULL; curr = curr->next){
        if (x <= _MAX && x >= _MAX - 100) {
            _MAX = _MAX + 1000;
            parents = (char *) realloc (parents, _MAX * sizeof(char));
        }
        sprintf(parent_string, "%ld:%ld,%ld:%ld", curr->left, curr->key, curr->right, curr->key);

        x = append_char_array(x, parent_string, parents);

        if (curr->next != NULL) parents[x++] = ',';
        else parents[x++] = '\n';
    }
    parents[x] = '\0';

    return x;
}


/*
 Find all of the subtrees in the tree
 */
void get_subtree_of_node(TreeNode* root, CharLink*& head) {
    char* buffer;
    // Check to see if it is an internal node (key < 0) or a leaf
    if (root->left == NULL && root->right == NULL) {
        buffer = (char*) malloc(100);
        for (int x = 0; x < 100; x++)
            buffer[x] = '\0';
        sprintf(buffer, "%ld", root->key);
        add_subtree(head, buffer);
        return;
    }
    get_subtree_of_node(root->right, head);
    CharLink* right_link = head;
    get_subtree_of_node(root->left, head);
    CharLink* left_link = head;

    // Join the last two subtrees
    int sum_length = get_char_array_length(right_link->subtree) + get_char_array_length(left_link->subtree) + 2;
    buffer = (char*) malloc(sum_length);
    for (int x = 0; x < sum_length; x++)
        buffer[x] = '\0';
    sprintf(buffer, "%s %s", right_link->subtree, left_link->subtree);
    add_subtree(head, buffer);
    return;
}

/*
 Find all of the subtrees in the tree
 */
void get_newick_subtrees(TreeNode* root, CharLink*& head) {
    char* buffer;
    // Check to see if it is an internal node (key < 0) or a leaf
    if (root->left == NULL && root->right == NULL) {
        buffer = (char*) malloc(100);
        for (int x = 0; x < 100; x++)
            buffer[x] = '\0';
        sprintf(buffer, "%ld", root->key);
        add_subtree(head, buffer);
        return;
    }
    get_newick_subtrees(root->right, head);
    CharLink* right_link = head;
    get_newick_subtrees(root->left, head);
    CharLink* left_link = head;

    // Join the last two subtrees
    int sum_length = get_char_array_length(right_link->subtree) + get_char_array_length(left_link->subtree) + 20;
    buffer = (char*) malloc(sum_length);
    for (int x = 0; x < sum_length; x++)
        buffer[x] = '\0';
    sprintf(buffer, "(%s,%s)%ld", right_link->subtree, left_link->subtree, root->key);
    add_subtree(head, buffer);
    return;
}


void get_subtree_of_node_helper(TreeNode* root, char *&subtrees, int &len_subtrees, const char delim) {
    CharLink * head = NULL;
    if (delim == ',')
        get_subtree_of_node(root, head);
    else
        get_newick_subtrees(root, head);
    int _MAX = 10000;
    subtrees = (char*) malloc(_MAX);
    // Parse the CharLink linked-list
    while (head) {
        if (len_subtrees <= _MAX && len_subtrees >= _MAX - 5000 ) {
            _MAX = _MAX + 10000;
            subtrees = (char *) realloc (subtrees, _MAX * sizeof(char));
        }
        len_subtrees = append_char_array(len_subtrees, head->subtree, subtrees);
        if (head->next) subtrees[len_subtrees++] = delim;
        head = head->next;
    }
    subtrees[len_subtrees] = '\0';
    deleteSubtreeList(head);
}


int get_node_relationships(char *tree_string, char *&children, char *&parents, char *&subtrees) {
    /*
    :param tree_string: tree_info['subtree_of_node'][node]
    Function loads the whole tree into a tree struct where each node has a child and a parent
    Then the tree is queried for its ONE PARENT and potentially MULTIPLE CHILDREN
    */

    // Step 1: Load the tree
    Link * linked_list = NULL;
    load_linked_list(tree_string, linked_list);
//    print_list(linked_list);

    // Step 2: Traverse the linked-list to get parents and children strings for each node
    int len_children = get_children_of_nodes(linked_list, children);
//    printf("Children:\n%s\n", children);
    int len_parents = get_parents_of_nodes(linked_list, parents);
//    printf("Parents:\n%s\n", parents);

    // Step 3: Convert the linked-list to a tree structure
    TreeNode* root = NULL;
    std::stack<TreeNode*> merge;
    load_tree_from_list(linked_list, root, merge);
    if (!merge.empty()) {
        std::cerr << "ERROR: Stack not empty after merging subtrees!" << std::endl;
        print_list(linked_list);
        while (!merge.empty()) {
            cout << "Not popped: " << merge.top()->key  << endl;
            merge.pop();
        }
        cout << tree_string << endl;
        return 0;
    }

//    print_tree(root);

    // Step 4: Traverse the tree to get all subtrees
    int len_subtrees = 0;
    get_subtree_of_node_helper(root, subtrees, len_subtrees, ',');

    //Step 5: Clean up the tree and linked list
    deleteTree(root);
    deleteList(linked_list);

    return len_children + len_parents + len_subtrees;
}


static PyObject *get_parents_and_children(PyObject *self, PyObject *args) {
    char* tree_string;
    if (!PyArg_ParseTuple(args, "s", &tree_string)) {
        return NULL;
    }

    char* children;
    char* parents;
    char* subtrees;

    int length = get_node_relationships(tree_string, children, parents, subtrees);
    if (length == 0)
        return Py_BuildValue("s", "$");

    children = (char *) realloc(children, (length + 3));

    int c_pos = 0;
    while (children[c_pos]){
        c_pos++;
    }

    int p_pos = 0;
    while (parents[p_pos])
        children[c_pos++] = parents[p_pos++];

    int t_pos = 0;
    while (subtrees[t_pos])
        children[c_pos++] = subtrees[t_pos++];
    children[c_pos] = '\0';

    free(parents);
    free(subtrees);

    return Py_BuildValue("s", children);
}


TreeNode* lca_helper(TreeNode* root, char** node_names, long& ancestor) {
    // TODO: Finish this function for potential use
    if (root == NULL) {
        return root;
    }
    int x = 0;
    bool all_contained = true;

    while (node_names[x]) {
        long query = atol(node_names[x]);
        cout << query << "\t";
        x++;
    }
    cout << endl;

    lca_helper(root->left, node_names, ancestor);
    lca_helper(root->right, node_names, ancestor);

    while (node_names[x]) {
        cout << root->key << endl;
        long query = atol(node_names[x]);
        if (!find(query, root))
            all_contained = false;
        x++;
    }
    cout << root->key << endl;

    if (all_contained && ancestor == 0) {
        ancestor = root->key;
        cout << root->key << endl;
    }
    return root;

}


static PyObject *lowest_common_ancestor(PyObject *self, PyObject *args) {
    char* tree_string;
    char* leaves_strung;
    long ancestor = 0;

    if (!PyArg_ParseTuple(args, "ss", &tree_string, &leaves_strung)) {
        return NULL;
    }

    // Get the node numbers to find LCA
    char **leaves = csv_to_list(leaves_strung);

    // Step 1: Load the tree
    Link * linked_list = NULL;
    load_linked_list(tree_string, linked_list);
    // Step 2: Convert the linked-list to a tree structure
    TreeNode* root = NULL;
    std::stack<TreeNode*> merge;
    load_tree_from_list(linked_list, root, merge);
    // Step 3: lca will return the root node of the LCA node for which all leaves are children
    lca_helper(root, leaves, ancestor);

    free(leaves);
    return Py_BuildValue("i", ancestor);

}


static PyObject *build_subtrees_newick(PyObject *self, PyObject *args) {
    /*
     Function to parse the rooted, assigned tree and find all subtrees of the inserted node
     Algorithm:
        1. Load the tree (load_linked_list and load_tree_from_list)
        2. Recursively build the subtrees from leaves to root in Newick format and load into list (get_newick_subtrees)
        3. Parse the linked list for each subtree, separating them by semicolons
        4. Return string to Python
     */
    char* tree_string;
    if (!PyArg_ParseTuple(args, "s", &tree_string)) {
        return NULL;
    }
    Link * linked_list = NULL;
    load_linked_list(tree_string, linked_list);

    TreeNode* root = NULL;
    std::stack<TreeNode*> merge;
    load_tree_from_list(linked_list, root, merge);
    if (!merge.empty()) {
        std::cerr << "ERROR: Stack not empty after merging subtrees!" << std::endl;
        print_list(linked_list);
        return 0;
    }
    char* subtrees;
    int len_subtrees = 0;

    get_subtree_of_node_helper(root, subtrees, len_subtrees, ';');

    deleteTree(root);
    deleteList(linked_list);

    return Py_BuildValue("s", subtrees);
}
