#ifndef SYMOP_H
#define SYMOP_H

#include <casm/external/Eigen/src/Core/Matrix.h>
#include <casmutils/xtal/site.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/sym/cartesian.hpp>
#include <string>
#include <vector>

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

struct VectorPeriodicCompare_f
{
       VectorPeriodicCompare_f(const Eigen::Vector3d& vector, double tol, const casmutils::xtal::Lattice& lattice);
       bool operator()(const Eigen::Vector3d& other) const;
    private:
        const Eigen::Vector3d m_vector;
        double m_precision;
        const casmutils::xtal::Lattice m_lattice;

};

struct SitePeriodicCompare_f
{
        SitePeriodicCompare_f(const casmutils::xtal::Site site, double tol, const casmutils::xtal::Lattice& lattice);
        bool operator()(casmutils::xtal::Site othersite) const;
  private:
        const casmutils::xtal::Site m_site;
        double m_precision;
        const casmutils::xtal::Lattice m_lattice;
};

struct CartOpCompare_f
{

    CartOpCompare_f(casmutils::sym::CartOp input1, double tol=1e-5);
    bool operator()(const casmutils::sym::CartOp& element2) const;

private:
    const casmutils::sym::CartOp element1;
    double tol;
};

//TODO: Better name. This doesn't describe anything except the function signature
class CartesianBinaryComparator_f
{
    public:
            CartesianBinaryComparator_f(double tol);
            bool operator()(const casmutils::sym::CartOp& element1, const casmutils::sym::CartOp& element2) const;
    private:
            double tol;
};

class BinaryCartOpPeriodicCompare_f
{
	public:
		BinaryCartOpPeriodicCompare_f(const casmutils::xtal::Lattice& lattice, double tol);
		bool operator()(const casmutils::sym::CartOp& element1, const casmutils::sym::CartOp& element2) const;
	private:
		double tol;
        casmutils::xtal::Lattice m_lattice;
};

class BinaryCartOpPeriodicMultiplier_f
{
    public:
        BinaryCartOpPeriodicMultiplier_f(const casmutils::xtal::Lattice& lattice, double tol);
        casmutils::sym::CartOp operator()(const casmutils::sym::CartOp& operation1, const casmutils::sym::CartOp& operation2) const;
    private:
        casmutils::xtal::Lattice m_lattice;
        double tol;
};

bool operator==(const casmutils::sym::CartOp& lhs, const casmutils::sym::CartOp& rhs);

#endif
