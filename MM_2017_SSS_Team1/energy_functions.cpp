#include <iostream>
#include <string>
#include <cmath>
#include <math.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>

namespace py = pybind11;

//double minimum_image_distance(py::array_t<double> r_i, py::array_t<double> r_j, double box_length)
//{
//  	py::buffer_info r_i_info = r_i.request();	
//	py::buffer_info r_j_info = r_j.request();
//        py::array_t<double> r_ix = r_i[0];
//        py::array_t<double> r_iy = r_i[1];
//        py::array_t<double> r_iz = r_i[2];
//        py::array_t<double> r_jx = r_j[0];
//        py::array_t<double> r_jy = r_j[1];
//        py::array_t<double> r_jz = r_j[2];
//              
//        double rij_x = r_ix - r_jx;
//        double rij_y = r_iy - r_jy;
//        double rij_z = r_iz - r_jz;
//        
//	r_2 = pow(rij_x,2.0) + pow(rij_y,2.0) + pow( rij_z,2.0);
//	r = pow(r_2, 0.5);
//
//	return r;
//}

//		py::array_t<double> rij = r_i - r_j;
//		rij = rij - box_length * int round(rij / box_length);
//		py::buffer_info rij_info = rij.request();
//		double dot = 0.0;
//		const double * r_ij_data = static_cast<double *>(r_ij_info.ptr);
//		for (size_t i = 0; i < r_ij_info.shape[0]; i++)
//		{
//			dot += r_ij_data[i] * r_ij_data[i];
// 		}
//		double rij = pow(dot,0.5)
//                return rij;

double lennard_jones_potential(double r)
{
    double sig_by_r6 = 1 / pow(r, 6.0);
    double sig_by_r12 = 1/ pow(sig_by_r6, 2.0);
    return 4.0 * (sig_by_r12 - sig_by_r6);
}
double lj_gradient_co(double r)
{
    double sig_by_r8 = 1 / pow(r, 8.0);
    double sig_by_r14 = 1 / pow(r,14.0);  
    return -48.0 * (sig_by_r14 - 0.5* sig_by_r8);
}

double switching_function(double r, double r_sw, double r_cut)
{
    	double r_2 = pow( r, 2.0);
	double r_sw2 = pow(r_sw, 2.0);
	double r_cut2 = pow(r_cut,2.0);
	if (r_2 <= r_sw2)
	{
		return 1.0;
	}
	else if (r_sw2 < r_2 and r_2 < r_cut2)
	{
		double s = pow(r_cut2 - r_2, 2.0) * (r_cut2 + 2 * r_2 - 3 * r_sw2) / pow(r_cut2 - r_sw2, 3.0);
		return s;
	}
	else
	{
		return 0.0;
	}
}

//sd::vector<double> deltar_SP(sd::vector<double> r_i, sd::vector<double> r_j)
//{
//	std::vector<double> rij = r_i - r_j
//	rij.size(); 
//	double r_2 = pow( r, 2.0);
//   	double r_sw2 = pow(r_sw, 2.0);
//    	double r_cut2 = pow(r_cut,2.0);		
//
//
//}




PYBIND11_PLUGIN(energy_functions)
{
    py::module m("energy_functions", "MM_1 basic module");
    m.def("lennard_jones_potential", &lennard_jones_potential);
    m.def("lj_gradient_co", &lj_gradient_co);
    m.def("switching_function", &switching_function);
    return m.ptr();
    
}
