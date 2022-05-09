/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "calibration.h"
#include "matrix_algo.h"

// added by the team
#include <iostream>
#include <fstream>


using namespace easy3d;



/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 * @return True on success, otherwise false. On success, the camera parameters are returned by
 */
bool Calibration::calibration(
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx, double& fy,    /// output: the focal length (in our slides, we use 'alpha' and 'beta'),
        double& cx, double& cy,    /// output: the principal point (in our slides, we use 'u0' and 'v0'),
        double& skew,              /// output: the skew factor ('-alpha * cot_theta')
        Matrix33& R,               /// output: the 3x3 rotation matrix encoding camera orientation.
        Vector3D& t)               /// output：a 3D vector encoding camera translation.
{

//    std::cout << "\n[Liangliang]:\n"
//                 "\tThe input parameters of this function are:\n"
//                 "\t\t- points_3d: An array of 3D points (input to this function)\n"
//                 "\t\t- points_2d: An array of 2D image points (input to this function)\n"
//                 "\tThis function must return either 'true' on success or 'false' otherwise. On success, the camera\n"
//                 "\tparameters are returned by the following variables:\n"
//                 "\t\t- fx and fy: the focal lengths (in our slides, we use 'alpha' and 'beta')\n"
//                 "\t\t- cx and cy: the principal point (in our slides, we use 'u0' and 'v0')\n"
//                 "\t\t- skew:      the skew factor ('-alpha * cot_theta')\n"
//                 "\t\t- R:         the 3x3 rotation matrix encoding camera orientation\n"
//                 "\t\t- t:         a 3D vector encoding camera location.\n"
//                 "\tIMPORTANT: don't forget to write your recovered parameters to the above variables." << std::endl;

    // TODO: check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    std::cout << "Checks for points >= 6" << std::endl;
    if (points_3d.size() >= 6 && points_2d.size() >= 6 && points_2d.size() == points_3d.size()){
        std::cout<< "the input is valid." << std::endl;
    }
    else{
        std::cout<< "the input is not valid." << std::endl;
        exit(1);
    }

    // TODO: construct the P matrix (so P * m = 0).
    // Define an m-by-12 double valued Projection matrix P.
    int m = 2 * points_2d.size(), n = 12;
    std::cout<<"m: "<<m <<std::endl;
    std::cout<<"n: "<<n <<std::endl;
    Matrix P(m, n, 0.0);
    int i0 = 0;
    int i = 0;
    int ii = 1;
    for (auto el3: points_3d){
        int j = 0;
        for (auto el2: points_2d){
            // I want something to be happening when i == j .
            if (i0 == j){
                // then I have everything I need to construct the row
                // set the first row by a 9-dimensional vector
                P.set_row(i, {el3.x(), el3.y(), el3.z(), 1, 0, 0, 0, 0, -el2.x()*el3.x(), -el2.x()*el3.y(), -el2.x()*el3.z(), -el2.x()}); // {....} is equivalent to a std::vector<double>
                P.set_row(ii, {0, 0, 0, 0, el3.x(), el3.y(), el3.z(), 1, -el2.y()*el3.x(), -el2.y()*el3.y(), -el2.y()*el3.z(), -el2.y()});
            }
            j++;
        }
        i= i+2;
        ii = ii+2;
        i0++;
    }
    //std::cout << "P: \n" << P << std::endl;


    // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.

    //Compute the SVD decomposition of A
    Matrix U(m, m, 0.0);
    Matrix S(m, 12, 0.0); // diagonal matrix that contains the singular values in descending order.
    Matrix V(12, 12, 0.0); // stores the corresponding singular vectors. We choose the singular vector p12 which minimises the error.
    svd_decompose(P, U, S, V);

    // PRINTS
//    std::cout << "U: \n" << U << std::endl;
//    std::cout << "S: \n" << S << std::endl;
//    std::cout << "V: \n" << V << std::endl;

    // The solution is the last column of V.
    // Create Matrix M.
    int r = 0;
    int c = 0;
    std::vector<double> em;
    std::cout << "the last column of V: " << std::endl;
    for(int k = 0; k < V.rows(); k++) {
        std::cout << V[k][V.rows()-1] << std::endl;
        em.emplace_back(V[k][V.rows()-1]);
        r++;
    }
    //std::cout << "elements for em: " << em << std::endl;

    // my Projection Matrix M
    Matrix34 M(em[0], em[1], em[2], em[3],
               em[4], em[5], em[6], em[7],
               em[8], em[9], em[10], em[11]);
    //std::cout << "M: \n" << M << std::endl;


    // the A(3x3) = the 3x3 elements of the M(4x4) Matrix
    Matrix33 A(em[0], em[1], em[2],
               em[4], em[5], em[6],
               em[8], em[9], em[10]);
    //std::cout << "A: \n" << A << std::endl;

    Matrix b(3, 1, 0.0);
    b[0][0] = em[3];
    b[1][0] = em[7];
    b[2][0] = em[11];
    //std::cout << "b: \n" << b << std::endl;


//    // TODO: extract intrinsic parameters from M.
    // Calculate the intrinsic parameters
    double p = 1/A.get_row(2).length(); // UNKNOWN SIGN --> NEED TO FIGURE OUT IF IT'S POSITIVE OR NEGATIVE
    double dota13 = dot(A.get_row(0),A.get_row(2));
    double dota23 = dot(A.get_row(1),A.get_row(2));
    Vector3D cra13 = cross(A.get_row(0),A.get_row(2));
    Vector3D cra23 = cross(A.get_row(1),A.get_row(2));
    cx = (p*p) * dota13 ;
    cy = (p*p) * dota23 ;
    double cos_theta = - dot(cra13,cra23) / (cra13.length() * cra23.length());   // HERE there might be a problem in case there is a zero denominator.
    double sin_theta = sqrt(1-cos_theta*cos_theta);
    fx = (p*p) * cra13.length() * sin_theta; //  fx = alpha in his notation
    fy = (p*p) * cra23.length(); // beta/sinΘ

    // PRINT to check the intrinsics
//    std::cout << "p: " << p << std::endl;
//    std::cout << "p^2: " << p*p << std::endl;
//    std::cout << "dota13: " << dota13 << std::endl;
//    std::cout << "dota23: " << dota23 << std::endl;
//    std::cout << "cra13: " << cra13 << std::endl;
//    std::cout << "cra23: " << cra23 << std::endl;
//    std::cout << "cx: " << cx << std::endl;
//    std::cout << "cy: " << cy << std::endl;
//    std::cout << "dot(cra13,cra23): " << dot(cra13,cra23) << std::endl;
//    std::cout << "cos_theta: " << cos_theta << std::endl;
//    std::cout << "sin_theta: " << sin_theta << std::endl;
//    std::cout << "fx: " << fx << std::endl;
//    std::cout << "fy: " << fy << std::endl;



//    // TODO: extract extrinsic parameters from M.

    Vector3D r1 = cra23 / cra23.length();
    Vector3D r3 = p * A.get_row(2);
    Vector3D r2 = cross(r3,r1);


    // define elements of Matrix K
    skew = - fx * cos_theta/sin_theta;

    //Matrix33 K(fx, skew, cx, 0, fy, cy, 0,0,1);
    Matrix33 K;
    K[0][0] = fx;
    K[0][1] = skew;
    K[0][2] = cx;
    K[1][0] = 0;
    K[1][1] = fy;
    K[1][2] = cy;
    K[2][0] = 0;
    K[2][1] = 0;
    K[2][2] = 1;
    std::cout << "K: \n" << K << std::endl;

//    // Compute the inverse of a matrix
    Matrix invK;
    inverse(K, invK);
    Matrix t1 = p * (invK * b);
    t[0] = t1[0][0];
    t[1] = t1[1][0];
    t[2] = t1[2][0];

    R.set_row(0, r1);
    R.set_row(1, r2);
    R.set_row(2, r3);

    std::cout << "R: " << R << std::endl;

    std::cout << "t: " << t << std::endl;


    // TODO by team: check if the Calibration is successful
    // Restructure the M(3x4) matrix with the help of K(3x3) and [R,T](3x4) matrices.
    Matrix34 Rt;
    Rt[0][0] = r1[0];
    Rt[0][1] = r1[1];
    Rt[0][2] = r1[2];
    Rt[0][3] = t[0];
    Rt[1][0] = r2[0];
    Rt[1][1] = r2[1];
    Rt[1][2] = r2[2];
    Rt[1][3] = t[1];
    Rt[2][0] = r3[0];
    Rt[2][1] = r3[1];
    Rt[2][2] = r3[2];
    Rt[2][3] = t[2];

    std::cout << "Rt: " << Rt << std::endl;

    // reconstruction of M matrix
    M = K * Rt;

    std::cout << "coordinates of 2D point: " << std::endl;
    for (const auto& p: points_3d){
        Vector4D q = p.homogeneous();
        //points_2d_new.emplace_back(M*p);
        Vector3D v = M * q; // M is 3 by 4
        Vector2D w = v.cartesian();
        std::cout << w << std::endl;
    }



//    std::cout << "\n\tTODO: After you implement this function, please return 'true' - this will trigger the viewer to\n"
//                 "\t\tupdate the rendering using your recovered camera parameters. This can help you to visually check\n"
//                 "\t\tif your calibration is successful or not.\n\n" << std::flush;
    return true;
} // end of bool function Calibration::calibration
