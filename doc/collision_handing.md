# Collision Handing

碰撞检测是一个从粗到细的过程。

## Broad-phase Collision Culling

-> BVH 对目标对象的划分

-> SP 对目标空间的划分

相比较下，SP更易于实现，且一些方法支持 GPU。

BVH的优点是当目标对象变化后，可以只需要对BV更新。

### Bounding Volume Hierarchies

粗筛选 (broad phase) 一般从构建包围盒开始（当然还有一些 2D 例如分离轴 // separating axis 方法等等）// bounding volumes

​	-> 外接圆/球 // circle and sphere

​	-> 轴对齐包围盒 // axis aligned bounding box

​	-> 定向包围盒 // oriented bounding box 

​	-> k-DOP，一般化的 AABB //k-sides discrete oriented polytope 

​	-> 凸包 // convex hull 

包围盒需要考虑的特性：

​	-> 紧密度 -> 加大碰撞可能性

​	-> 构建方法 -> 算法效率

​	-> 相交测试 -> 算法效率

​	-> 全局变换

​	-> 空间占用

为了进一步得到准确的碰撞结果，需要结构化输入以加快区域搜索。

-> 树状结构化划分是常见的做法  // hierarchical clustering strategy

​	-> 自顶向下 // top-down

​		-> 分离轴 // splitting axis -> 最长轴等

​		-> 分离点 // splitting point -> 平均值，中值等

​	-> 自底向上 // bottom-up 

​		-> 子节点的选取 //leaf object clusters

​		-> 临近点计算 // merge children into parent

### Spatial Partitionings / Spatial Hashing 

方法分类：

​	-> 归一化细分 square or cube //uniform subdivision->将空间分隔成相同的正方形 (2D) 或正方体 (3D)

​	-> 四叉树或八叉树 // quadtree for 2d, octree for 3D 

​	-> kD-tree 

​	-> 二叉空间分割 // binary space partitioning

空间分割方法的痛点在于当目标对象处于分割线中间时的处理。

## Narrow-Phase Collision Test

通过粗筛选得到或是BV或SH结构的碰撞对，需要进一步对结构里的对象做测试。

不同的处理方式分类：

​	-> 连续碰撞检测 // continue collision detection -> 多用于模拟仿真 -> 考虑物理量下的碰撞检测 -> 加入到最小化能量

​	-> 离散碰撞检测 // discrete collision detection -> 静态的碰撞检测 -> 考虑各种复形的相交情况

### Discrete Collision Detection

