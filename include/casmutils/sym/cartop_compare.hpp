#ifndef SYMOP_H
#define SYMOP_H

#include "../../submodules/eigen-git-mirror/Eigen/Core"
#include "../../submodules/eigen-git-mirror/Eigen/Dense"
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/sym/cartesian.hpp>
#include <string>
#include <vector>
#include "lattice.hpp"

//class CartOp
//{
//public:
//    CartOp(const Eigen::Matrix3d& cart_matrix, const Eigen::Vector3d& translation={0, 0, 0});
//    Eigen::Vector3d get_translation() const;
//    Eigen::Matrix3d get_cart_matrix() const;
//
//    Eigen::Matrix3d m_cart_matrix;
//private:
//
//    Eigen::Vector3d m_translation;
//};
//
//CartOp operator*(const CartOp& lhs, const CartOp& rhs);


//TODO: Refactor. Just elimitate this, we only want the binary comparisons
//If you want unary comparisons, make a fancy template
//Change later
struct CartOpCompare_f
{

    CartOpCompare_f(CartOp input1, double tol=1e-5);
    bool operator()(const CartOp& element2) const;

private:
    const CartOp element1;
    double tol;
};

//TODO: Better name. This doesn't describe anything except the function signature
class CartesianBinaryComparator_f
{
    public:
            CartesianBinaryComparator_f(double tol);
            bool operator()(const CartOp& element1, const CartOp& element2) const;
    private:
            double tol;
};

class BinaryCartOpPeriodicCompare_f
{
	public:
		BinaryCartOpPeriodicCompare_f(const Lattice& lattice, double tol);
		bool operator()(const CartOp& element1, const CartOp& element2) const;
	private:
		double tol;
		Lattice m_lattice;
};

class BinaryCartOpPeriodicMultiplier_f
{
    public:
        BinaryCartOpPeriodicMultiplier_f(const Lattice& lattice, double tol);
        CartOp operator()(const CartOp& operation1, const CartOp& operation2) const;
    private:
        Lattice m_lattice;
        double tol;
};

bool operator==(const CartOp& lhs, const CartOp& rhs);

#endif
