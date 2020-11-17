#include <casmutils/sym/cartesian.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/site.hpp>
#include <casmutils/xtal/lattice.hpp>

//
//CartOp operator*(const CartOp& lhs, const CartOp& rhs)
//{
//    Eigen::Vector3d translation = lhs.matrix * rhs.translation + lhs.translation;
//    Eigen::Matrix3d product = lhs.matrix * rhs.matrix;
//    CartOp symop_product(product, translation);
//    return symop_product;
//}
//
//
//TODO: Make this disappear
CartOpCompare_f::CartOpCompare_f(casmutils::sym::CartOp input1, double tol) : element1(input1), tol(tol) {}

bool CartOpCompare_f::operator()(const casmutils::sym::CartOp& element2) const
{
    return (element1.matrix.isApprox(element2.matrix, tol) &&
            (element1.translation.isApprox(element2.translation, tol)));
}

CartesianBinaryComparator_f::CartesianBinaryComparator_f(double tol) : tol(tol) {}

bool CartesianBinaryComparator_f::operator()(const casmutils::sym::CartOp& element1, const casmutils::sym::CartOp& element2) const
{

   return (element1.matrix.isApprox(element2.matrix, tol) &&
            (element1.translation.isApprox(element2.translation, tol)));
}


BinaryCartOpPeriodicCompare_f::BinaryCartOpPeriodicCompare_f(const casmutils::xtal::Lattice& lattice, double tol) : m_lattice(lattice), tol(tol) {}

bool BinaryCartOpPeriodicCompare_f::operator()(const casmutils::sym::CartOp& element1, const casmutils::sym::CartOp& element2) const
{    
    casmutils::xtal::Site temp_site1 = casmutils::xtal::Site( casmutils::xtal::Coordinate(element1.translation),std::string("xx") );
    casmutils::xtal::Site temp_site2 = casmutils::xtal::Site(casmutils::xtal::Coordinate(element2.translation),std::string("xx"));
    SitePeriodicCompare_f translation_comparison(temp_site1, tol, m_lattice);
    
    casmutils::sym::CartOp symop1(element1.matrix);
    casmutils::sym::CartOp symop2(element2.matrix);
	CartesianBinaryComparator_f compare(tol);

	return compare(symop1, symop2) && translation_comparison(temp_site2);
}

BinaryCartOpPeriodicMultiplier_f::BinaryCartOpPeriodicMultiplier_f(const Lattice& lattice, double tol) : m_lattice(lattice), tol(tol) {}

casmutils::sym::CartOp BinaryCartOpPeriodicMultiplier_f::operator()(const casmutils::sym::CartOp& operation1, const casmutils::sym::CartOp& operation2) const 
{
    casmutils::sym::CartOp full_operation_product = operation1 * operation2;
    Eigen::Vector3d op_product_periodic_tranlation = bring_within(m_lattice, tol, full_operation_product.translation);
    casmutils::sym::CartOp final_product(full_operation_product.matrix, op_product_periodic_tranlation);
    return final_product;
}

bool operator==(const casmutils::sym::CartOp& lhs, const casmutils::sym::CartOp& rhs)
{
    CartesianBinaryComparator_f binarycompare(1e-6);
    return binarycompare(lhs, rhs);
}

