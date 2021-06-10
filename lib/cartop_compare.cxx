#include <casmutils/sym/cartesian.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/site.hpp>
#include <casmutils/xtal/lattice.hpp>
#include "../include/casmutils/sym/cartop_compare.hpp"

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

//Eigen::Vector3d bring_within(casmutils::xtal::Lattice lattice, double tol, Eigen::Vector3d vector){
//    return vector;
//}

VectorPeriodicCompare_f::VectorPeriodicCompare_f(const Eigen::Vector3d& vector, double tol, const casmutils::xtal::Lattice& lattice) : m_vector(vector), m_precision(tol), m_lattice(lattice){}

bool VectorPeriodicCompare_f::operator()(const Eigen::Vector3d& other) const
{
    Eigen::Vector3d vector1 = casmutils::xtal::cartesian_to_fractional( m_vector, m_lattice);
    Eigen::Vector3d vector2 = casmutils::xtal::cartesian_to_fractional(other, m_lattice);
    Eigen::Vector3d distance_vector = vector1-vector2;

    for (int i=0; i<distance_vector.size(); i++){
        distance_vector(i) = distance_vector(i) - std::round(distance_vector(i));
    }
    Eigen::Vector3d cartesian_distance_vector = m_lattice.column_vector_matrix()*distance_vector;
    return std::abs(cartesian_distance_vector.norm())<m_precision;
}

SitePeriodicCompare_f::SitePeriodicCompare_f(const casmutils::xtal::Site site, double tol, const casmutils::xtal::Lattice& lattice) : m_site(site), m_precision(tol), m_lattice(lattice){
}
bool SitePeriodicCompare_f::operator()(casmutils::xtal::Site othersite) const
{
    VectorPeriodicCompare_f compare(m_site.cart(), m_precision, m_lattice);
    return m_site.label()==othersite.label() &&compare(othersite.cart());
}

CartOpCompare_f::CartOpCompare_f(casmutils::sym::CartOp input1, double tol) : element1(input1), tol(tol) {}

bool CartOpCompare_f::operator()(const casmutils::sym::CartOp& element2) const
{
    return ((element1.matrix.isApprox(element2.matrix, tol) &&
            element1.translation.isApprox(element2.translation, tol)) &&
            (element1.is_time_reversal_active==element2.is_time_reversal_active));
}

CartesianBinaryComparator_f::CartesianBinaryComparator_f(double tol) : tol(tol) {}

bool CartesianBinaryComparator_f::operator()(const casmutils::sym::CartOp& element1, const casmutils::sym::CartOp& element2) const
{

   return (element1.matrix.isApprox(element2.matrix, tol) &&
            (element1.translation.isApprox(element2.translation, tol)) && 
            (element1.is_time_reversal_active==element2.is_time_reversal_active));
}


BinaryCartOpPeriodicCompare_f::BinaryCartOpPeriodicCompare_f(const casmutils::xtal::Lattice& lattice, double tol) : m_lattice(lattice), tol(tol) {}

bool BinaryCartOpPeriodicCompare_f::operator()(const casmutils::sym::CartOp& element1, const casmutils::sym::CartOp& element2) const
{    
    casmutils::xtal::Site temp_site1 = casmutils::xtal::Site( element1.translation, std::string("xx") );
    casmutils::xtal::Site temp_site2 = casmutils::xtal::Site(element2.translation, std::string("xx"));
    SitePeriodicCompare_f translation_comparison(temp_site1, tol, m_lattice);
    
    casmutils::sym::CartOp symop1(element1.matrix,{0,0,0},false);
    casmutils::sym::CartOp symop2(element2.matrix,{0,0,0}, false);
	CartesianBinaryComparator_f compare(tol);

	return compare(symop1, symop2) && translation_comparison(temp_site2) && (element1.is_time_reversal_active==element2.is_time_reversal_active);
}

BinaryCartOpPeriodicMultiplier_f::BinaryCartOpPeriodicMultiplier_f(const casmutils::xtal::Lattice& lattice, double tol) : m_lattice(lattice), tol(tol) {}

casmutils::sym::CartOp BinaryCartOpPeriodicMultiplier_f::operator()(const casmutils::sym::CartOp& operation1, const casmutils::sym::CartOp& operation2) const 
{
    casmutils::sym::CartOp full_operation_product = operation1 * operation2;
    Eigen::Vector3d op_product_periodic_tranlation = casmutils::xtal::bring_within_lattice(casmutils::xtal::cartesian_to_fractional(full_operation_product.translation, m_lattice), m_lattice);
    casmutils::sym::CartOp final_product(full_operation_product.matrix, op_product_periodic_tranlation, (operation1.is_time_reversal_active != operation2.is_time_reversal_active));
    return final_product;
}

bool operator==(const casmutils::sym::CartOp& lhs, const casmutils::sym::CartOp& rhs)
{
    CartesianBinaryComparator_f binarycompare(1e-6);
    return binarycompare(lhs, rhs);
}

