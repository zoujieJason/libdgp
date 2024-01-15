# KDTree

## Binary Search Tree

### How to build one?

 给定一组一维数据 -> 根据二分排序插入每个数据

```c++
void insert(Tree *tree, node *node)
{
  if reach a leaf of tree 
    node set to be the leaf and return 
  tree->value < node->value ? insert(tree->leftChild, node): insert(tree->rightChild, node);
}

Tree* build(List list)
{
  fvalue = the first value of the list;
  build a tree with fvalue;
	for value in list
    build a node with value;
    insert node into tree; 
  return tree;
}

```

->这种简单的插入式构建在极端情况下会退化树为链表 -> 完全取决于数据的输入序列

#### How to improve it?

-> 对输入序列排序 -> $O(n\log n)$ -> 后续的插入和删除操作无法保证平衡树

-> Self-balancing BST (AVL tree) -> 每个节点存储了平衡因子 // balance factor -> 需满足条件：
$$
\forall \text{node}\in \text{Tree}, \text{bf}(\text{node})=\text{Height}(\text{node}.\text{leftChild})-\text{Height}(\text{node}.\text{rightChild})\\
\text{s.t}\quad \text{bf}(\text{node})\in \{-1, 0, 1\}
$$
-> 树的插入、删除操作都需要满足上述条件

这里不再进一步叙述 AVL tree 的一些属性和实现。

当然还有一些其他的方法使树保持平衡的方法 //rebalance 

### Traverse

Search the value:

```c++

Value min_diff 
void traverse(Tree *node, Value value)
{
	if node is null
    return 
  min_diff = MIN(abs(node->value - value), min_diff)
  node->value < value ? traverse(node->leftChild, value): 
  traverse(node->rightChild, value);
}
```

average case: $O(\log n)$ -> worst case: $O(n)$ 

对一维的 BST 每一次迭代都能将值划分至所属的区间，所以无需额外的回溯来重新计算“最近值”。

## How it works?

BST的节点 (node) 一般既是**划分节点**，又是**数据节点**（当然可以分离设计）-> 处理多维，树节点只能存储*划分信息*，所有数据存储在叶子节点里

### Build

How to build a kdtree ? -> 给一组三维数据 -> 选择合适的超平面，分割左右子树数据 -> 迭代循环

```c++
void build(Tree *tree, List list)
{
  if the length of the list less than the threshold
    set as a leaf and return 
  find a hyperplane to split the list into part1 and part2;
  build a node of the tree with the hyperplane;
  build(tree->leftChild, list_part1);
  build(tree->rightChild, list_part2);
}
```

通常采用中位数分割方法 (median)，所以整个构建的复杂度为 $O(n\log n)$。

建树过程跟二叉树别无二致。

### Traverse

```c++
void depth_first_search(Node *node)
{
  if(node is nullptr)
    return
  depth_first_search(node->left); 
  depth_first_search(node->right);
}
```

整个遍历的过程，是深度优先寻找一个初值，然后在返回的过程中剪枝。

```c++
min_sq_dis; 
void traverse(Node * node, Point q)
{
  if(node is leaf)
  {
    min_sq_dis = distance(node, q); 
  }
  else 
  {
    if(split(node->hyperplane, q))
    {
      traverse(node->leftChild, q); 
      if is possible to find the closest point
      {
        traverse(node->rightChild, q); 
      }
    }
    else 
    {
      traverse(node->rightChild, q);
      if is possible to find the closest point
      {
        traverse(node->leftChild, q); 
      }
    }
  }
}
```

### Some Experiments

经过一番设计，更精细的剪枝计算会带来更大的计算代价，反倒不如一些简单的剪枝，如超平面二分。
