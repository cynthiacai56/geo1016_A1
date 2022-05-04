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
//    // TO BE REMOVED - BEGINNING
//    std::cout << "\nTODO: I am going to implement the calibration() function in the following file:\n"
//                 "\t    - calibration_method.cpp\n\n";
//
//    std::cout << "[Liangliang]:\n"
//                 "\tCamera calibration requires computing the SVD and inverse of matrices.\n"
//                 "\tIn this assignment, I provide you with a 'Matrix' and a 'Vector' data structures for storing and\n"
//                 "\tmanipulating matrices and vectors of arbitrary sizes. I also wrote some code to show you how to:\n"
//                 "\t    - compute the SVD of a matrix;\n"
//                 "\t    - compute the inverse of a matrix;\n"
//                 "\t    - compute the transpose of a matrix.\n\n"
//                 "\tFeel free to use any of the provided data structures and functions. The commonly used linear algebra\n"
//                 "\tfunctions are provided in the following files:\n"
//                 "\t    - Calibration/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
//                 "\t    - Calibration/vector.h  Vectors of arbitrary dimensions and related functions.\n"
//                 "\t    - Calibration/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
//                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
//                 "\tIn your final submission, please\n"
//                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
//                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
//                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
//                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n" << std::flush;
//
//    /// Below are a few examples showing some useful data structures and functions.
//
//    // This is a 1D array of 'double' values. Alternatively, you can use 'double mat[25]' but you cannot change it
//    // length. With 'std::vector', you can append/delete/insert elements, and much more. The 'std::vector' can store
//    // not only 'double', but also any other types of objects. In case you may want to learn more about 'std::vector'
//    // check here: https://en.cppreference.com/w/cpp/container/vector
//    std::vector<double> array = {1, 3, 3, 4, 7, 6, 2, 8, 2, 8, 3, 2, 4, 9, 1, 7, 3, 23, 2, 3, 5, 2, 1, 5, 8, 9, 22};
//    array.push_back(5); // append 5 to the array (so the size will increase by 1).
//    array.insert(array.end(), 10, 3);  // append ten 3 (so the size will grow by 10).
//
//    /// To access the value of an element.
//    double a = array[2];
//
//    /// define a 2D vector/point
//    Vector2D b(1.1, 2.2);
//
//    /// define a 3D vector/point
//    Vector3D c(1.1, 2.2, 3.3);
//
//    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
//    Vector2D p = c.cartesian();
//
//    /// get the Homogeneous coordinates of p
//    Vector3D q = p.homogeneous();
//
//    /// the length of a vector
//    double len = p.length();
//    /// the squared length of a vector
//    double sqr_len = p.length2();
//
//    /// the dot product of two vectors
//    double dot_prod = dot(p, q);
//
//    /// the cross product of two vectors
//    Vector cross_prod = cross(c, q);
//
//    /// normalize this vector
//    cross_prod.normalize();
//
//    // Define an m-by-n double valued matrix.
//    // Here I use the above array to initialize it. You can also use A(i, j) to initialize/modify/access its elements.
//    const int m = 6, n = 5;
//    Matrix A(m, n, array.data());    // 'array.data()' returns a pointer to the array.
////    std::cout << "M: \n" << A << std::endl;
//
//    /// define a 3 by 4 matrix (and all elements initialized to 0.0)
//    Matrix M(3, 4, 0.0);
//
//    /// set first row by a vector
//    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));
//
//    /// set second column by a vector
//    M.set_column(1, Vector3D(5.5, 5.5, 5.5));
//
//    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
//    Matrix33 B;
//
//    /// define and initialize a 3 by 3 matrix
//    Matrix33 T(1.1, 2.2, 3.3,
//               0, 2.2, 3.3,
//               0, 0, 1);
//
//    /// define and initialize a 3 by 4 matrix
//    Matrix34 P(1.1, 2.2, 3.3, 0,
//               0, 2.2, 3.3, 1,
//               0, 0, 1, 1);
//
//    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
//    Matrix W(15, 9, 0.0);
//    /// set the first row by a 9-dimensional vector
//    W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>
//
//    /// get the number of rows.
//    int num_rows = W.rows();
//
//    /// get the number of columns.
//    int num_cols = W.cols();
//
//    /// get the element at row 1 and column 2
//    double value = W(1, 2);
//
//    /// get the last column of a matrix
//    Vector last_column = W.get_column(W.cols() - 1);
//
//    /// define a 3 by 3 identity matrix
//    Matrix33 I = Matrix::identity(3, 3, 1.0);
//
//    /// matrix-vector product
//    Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4
//
//    Matrix U(m, m, 0.0);   // initialized with 0s
//    Matrix S(m, n, 0.0);   // initialized with 0s
//    Matrix V(n, n, 0.0);   // initialized with 0s
//
//    // Compute the SVD decomposition of A
//    svd_decompose(A, U, S, V);
//
//    // Now let's check if the SVD result is correct
//
//    // Check 1: U is orthogonal, so U * U^T must be identity
////    std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;
//
//    // Check 2: V is orthogonal, so V * V^T must be identity
////    std::cout << "V*V^T: \n" << V * transpose(V) << std::endl;
//
//    // Check 3: S must be a diagonal matrix
////    std::cout << "S: \n" << S << std::endl;
//
//    // Check 4: according to the definition, A = U * S * V^T
////    std::cout << "M - U * S * V^T: \n" << A - U * S * transpose(V) << std::endl;
//
//    // Compute the inverse of a matrix
//    Matrix invT;
//    inverse(T, invT);
//    // Let's check if the inverse is correct
////    std::cout << "B * invB: \n" << B * invB << std::endl;
//
//    // TODO: the above code just demonstrates some useful data structures and APIs. Please remove all above code in your
//    //       final submission.
//
//    //--------------------------------------------------------------------------------------------------------------
    // TO BE REMOVED - END


    // implementation starts ...

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
    // Define an m-by-n double valued Projection matrix P.
    const int m = 2 * points_2d.size(), n = 12;
    Matrix P(m, n, 0.0);

    // Now we have to set the elements equal to meaningful values.
    // Have a loop for every point stored in the points_2d and points_3d vectors.
    // and construct the 2 rows for the x, and y coordinate respectively.
    // construct the row, then push it to the matrix.
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
                P.set_row(i, {-el3.x(), -el3.y(), -el3.z(), -1, 0, 0, 0, 0, el2.x()*el3.x(), el2.x()*el3.y(), el2.x()*el3.z(), el2.x()}); // {....} is equivalent to a std::vector<double>
                P.set_row(ii, {0, 0, 0, 0, -el3.x(), -el3.y(), -el3.z(), -1, el2.y()*el3.x(), el2.y()*el3.y(), el2.y()*el3.z(), el2.y()});
            }
            j++;
        }
        i= i+2;
        ii = ii+2;
        i0++;
    }

    std::cout << "P: \n" << P << std::endl;


    // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.

    //Compute the SVD decomposition of A
    Matrix U(m, 12, 0.0);
    Matrix S(12, 12, 0.0); // diagonal matrix that contains the singular values in descending order.
    Matrix V(12, 12, 0.0); // stores the corresponding singular vectors. We choose the singular vector p12 which minimises the error.
    svd_decompose(P, U, S, V);

    std::cout << "U: \n" << U << std::endl;
    std::cout << "S: \n" << S << std::endl;
    std::cout << "V: \n" << V << std::endl;

    // The solution is the last column of VT. Get the last column of VT (transpose(V)). Which is the last row of V.
    // according to the definition, A = U * S * V^T
    // get the last row of V.
    // ... add them as elements to Projection Matrix M.
    // Create Matrix M.
    int r = 0;
    int c = 0;
    std::vector<double> em;
    for(int k = 0; k < V.rows(); k++) {
        //take the last row
        //std::cout << "the last row of V: " << V[V.rows()-1][k] << std::endl; //Last row
        //em.emplace_back(V[V.rows()-1][k]);
        // take the last column
        std::cout << "the last row of V: " << V[k][V.rows()-1] << std::endl; //Last row
        em.emplace_back(V[k][V.rows()-1]);
        r++;
    }
    std::cout << "elements for M: " << em << std::endl;

    // my Projection Matrix M
    Matrix34 M(em[0], em[1], em[2], em[3],
               em[4], em[5], em[6], em[7],
               em[8], em[9], em[10], em[11]);
    std::cout << "M: \n" << M << std::endl;

    // the H(3x3) = KR (the 3x3 elements of the M(4x4) Matrix).
    Matrix33 H(em[0], em[1], em[2],
               em[4], em[5], em[6],
               em[8], em[9], em[10]);
    //std::cout << "H: \n" << H << std::endl;

    Matrix hs(3, 1, 0.0);
    hs[0][0] = em[3];
    hs[1][0] = em[7];
    hs[2][0] = em[11];
    //std::cout << "h: \n" << hs << std::endl;

    // inverse of H
    Matrix invH = inverse(H);
    // Let's check if the inverse is correct
    std::cout << "H * invH: \n" << H * invH << std::endl; // Correct Inversion --> leads to I matrix.

    // Matrix Xo; // remember it must be homogeneous. so I need to go to Euclidean.
    /// matrix-matrix multiplication
    Matrix operator*(const Matrix &, const Matrix &);
    Matrix Xo = -operator*(invH, hs);
    //std::cout << "Xo: \n" << Xo << std::endl;

    // QR decomposition of matrix H into K(orthogonal matrix) and R(upper triangular matrix). H = K*R.
    // We follow the Gram-Schmidt process to implement the QR decomposition of matrix H.
    // STEPS: we retrieve the columns of matrix K. We save the retrieved column into a variable u.
    // based on variable u we define the k element of K matrix.
    /// STEP: traverse the columns of a matrix
    // EXAMPLE MATRIX - NOT FOR USE
//    Matrix H1(3, 3, 0.0);
//    H1[0][0] = 1;  H1[1][0] = 1;  H1[2][0] = 0;// H1[3][0] = 1;
//    H1[0][1] = 1;  H1[1][1] = 0;  H1[2][1] = 1; //H1[3][1] = -1;
//    H1[0][2] = 0;  H1[1][2] = 1;  H1[2][2] = 1; //H1[3][2] = 0;
//    std::cout << "H1: \n" << H1 << std::endl;

    Vector a1 = H.get_column(0);
    Vector e1 = a1/a1.length();
    Matrix33 K;
    K.set_column(0, e1);

    for (int cl = 1; cl < H.cols(); cl++){
        Vector acur = H.get_column(cl);
        Vector3D subtract(0,0,0);
        for (auto index=0; index < K.cols(); index++){
            subtract = subtract + (K.get_column(index))*(dot(acur, (K.get_column(index))));
        }
        Vector ucur =  acur - subtract;
        Vector ecur = ucur/ucur.length();
        K.set_column(cl, ecur);
    }
    //std::cout << "K: \n" << K << std::endl;

    R[0][0] = dot(H.get_column(0), K.get_column(0));
    R[0][1] = dot(H.get_column(1), K.get_column(0));
    R[0][2] = dot(H.get_column(2), K.get_column(0));
    R[1][0] = 0;
    R[1][1] = dot(H.get_column(1), K.get_column(1));
    R[1][2] = dot(H.get_column(2), K.get_column(1));
    R[2][0] = 0;
    R[2][1] = 0;
    R[2][2] = dot(H.get_column(2), K.get_column(2));
    //std::cout << "R: \n" << R << std::endl;

    // FLIP the matrices:
    Matrix33 FlipRzp(-1,0,0,0,-1,0,0,0,1);
    K = operator*(K,FlipRzp);
    R = operator*(FlipRzp, R);
    //std::cout << "K new: \n" << K << std::endl;
    std::cout << "R new: \n" << R << std::endl;


    // fx, fy,    /// output: the focal length (in our slides, we use 'alpha' and 'beta'),
    // cx, cy,    /// output: the principal point (in our slides, we use 'u0' and 'v0'),
    // skew,      /// output: the skew factor ('-alpha * cot_theta')
    // R,               /// output: the 3x3 rotation matrix encoding camera orientation.
    // Vector3D& t      /// output：a 3D vector encoding camera translation.

    fx = K[0][0];
    //std::cout << "fx: " << fx << std::endl;
    cx = K[0][2];
    //std::cout << "cx: " << cx << std::endl;
    cy = K[1][2];
    //std::cout << "cy: " << cy << std::endl;
    fy = K[1][1]; // ARE WE USING THE SIMPLE VERSION OF K MATRIX?
    //std::cout << "fy: " << fy << std::endl;
    skew = K[0][1];
    //std::cout << "skew: " << skew << std::endl;
    t = Vector3D(*Xo[0], *Xo[1], *Xo[2]);
    //std::cout << "t: " << t << std::endl;



    // TODO: extract intrinsic parameters from M.

    // TODO: extract extrinsic parameters from M.


    // TODO by team: check if the Calibration is successful
    // points_3d
    //std::vector<Vector2D> points_2d_new;
//    std::cout << "coordinates of 2D point: " << std::endl;
//    for (const auto& p: points_3d){
//        Vector4D q = p.homogeneous();
//        //points_2d_new.emplace_back(M*p);
//        Vector3D v = M * q; // M is 3 by 4
//        Vector2D w = v.cartesian();
//        std::cout << w << std::endl;
//    }

    Matrix34 Rt;
    Rt[0][0] = R[0][0];
    Rt[0][1] = R[0][1];
    Rt[0][2] = R[0][2];
    Rt[0][3] = *Xo[0];
    Rt[1][0] = R[1][0];
    Rt[1][1] = R[1][1];
    Rt[1][2] = R[1][2];
    Rt[1][3] = *Xo[1];
    Rt[2][0] = R[2][0];
    Rt[2][1] = R[2][1];
    Rt[2][2] = R[2][2];
    Rt[2][3] = *Xo[2];
    std::cout << "Rt: \n" << Rt << std::endl;
//
    for (const auto& p: points_3d){
        Vector4D q = p.homogeneous();
        //points_2d_new.emplace_back(M*p);
        Vector3D v = (operator*(K,Rt)) * q; // M is 3 by 4
        Vector2D w = v.cartesian();
        std::cout << w << std::endl;
    }

//    std::cout << "\n\tTODO: After you implement this function, please return 'true' - this will trigger the viewer to\n"
//                 "\t\tupdate the rendering using your recovered camera parameters. This can help you to visually check\n"
//                 "\t\tif your calibration is successful or not.\n\n" << std::flush;
    return false;
} // end of bool function Calibration::calibration



















