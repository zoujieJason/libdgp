#include <dgp/kdtree.h>

#include <iostream>

#include <chrono>

int main(int argc, char **argv)
{
    const int dim = std::stoi(argv[1]); 
    const int n_samplings = std::stoi(argv[2]);
    const int n_random = std::stoi(argv[3]);

    std::srand((unsigned int) time(0));
    Eigen::initParallel();
    using Scalar = float;
    using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>; 
    using VectorType = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;

    MatrixType V = MatrixType::Random(n_samplings, dim);
    dgp::KDTree tree(V);
    for(int i = 0; i < n_random; ++i)
    {
        const VectorType q = VectorType::Random(dim);
        Scalar sq; 
        VectorType p;
        std::cout << i << " ";
        {
            auto start = std::chrono::high_resolution_clock::now();
            
            tree.brutal_search(q, p, sq);

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<Scalar> elapsed = end - start;
            std::cout << elapsed.count() << " ";
        }

        int cnt0; 
        {
            Scalar dist_sq; 
            VectorType closet_p;
            auto start = std::chrono::high_resolution_clock::now();

            cnt0 = tree.recursive_search(q, closet_p, dist_sq);

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<Scalar> elapsed = end - start;
            std::cout << elapsed.count() << " ";
            //std::cout << sq-dist_sq << " " << (p-closet_p).squaredNorm() << " " << cnt0 << " ";
        }

        int cnt1; 
        {
            Scalar dist_sq; 
            VectorType closet_p;
            auto start = std::chrono::high_resolution_clock::now();

            cnt1 = tree.iterative_search(q, closet_p, dist_sq);

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<Scalar> elapsed = end - start;
            std::cout << elapsed.count() << " ";
            //std::cout << sq-dist_sq << " " << (p-closet_p).squaredNorm() << " " << cnt1 << " ";
        }

        int cnt2; 
        {
            Scalar dist_sq; 
            VectorType closet_p;
            auto start = std::chrono::high_resolution_clock::now();

            cnt2 = tree.recursive_pruning_search(q, closet_p, dist_sq);

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<Scalar> elapsed = end - start;
            std::cout << elapsed.count() << " ";
            //std::cout << sq-dist_sq << " " << (p-closet_p).squaredNorm() << " " << cnt2 << " ";
        }
        std::cout << cnt0 - cnt1 <<  " " << static_cast<Scalar>(cnt0 - cnt2) / static_cast<Scalar>(cnt0) << std::endl;
    }
    return 0;
}