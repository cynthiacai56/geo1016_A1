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

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>

// -- added by the team --
#include <random>


using namespace easy3d;

// declare the functions that follow the triangulation function -- created by the team
Matrix33 transf_matrix(const std::vector<Vector2D>& im_points);
Matrix34 construct_Rt(const Matrix33& R, const Vector3D& t);
int count_positive_z(const std::vector<Vector2D>& img0, const std::vector<Vector2D>& img1, const Matrix34& KRt);
std::vector<Vector3D> reconstruct3Dpoints(const std::vector<Vector2D>& img0, const std::vector<Vector2D>& img1, const Matrix34& KRt, const Matrix34& M0);
int count_positive_z(const std::vector<Vector3D>& points3d);



/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */
bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for the sub-tasks. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or put them in one or multiple separate files.

    std::cout << "\nTODO: I am going to implement the triangulation() function in the following file:" << std::endl
              << "\t    - triangulation_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tFeel free to use any provided data structures and functions. For your convenience, the\n"
                 "\tfollowing three files implement basic linear algebra data structures and operations:\n"
                 "\t    - Triangulation/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/vector.h  Vectors of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please\n"
                 "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations.\n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n" << std::flush;

//    /// Below are a few examples showing some useful data structures and APIs.
//
//    /// define a 2D vector/point
//    Vector2D b(1.1, 2.2);
//
//    /// define a 3D vector/point
//    Vector3D a(1.1, 2.2, 3.3);
//
//    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
//    Vector2D p = a.cartesian();
//
//    /// get the Homogeneous coordinates of p
//    Vector3D q = p.homogeneous();
//
//    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
//    Matrix33 A;
//
//    /// define and initialize a 3 by 3 matrix
//    Matrix33 T(1.1, 2.2, 3.3,
//               0, 2.2, 3.3,
//               0, 0, 1);
//
//    /// define and initialize a 3 by 4 matrix
//    Matrix34 M(1.1, 2.2, 3.3, 0,
//               0, 2.2, 3.3, 1,
//               0, 0, 1, 1);
//
//    /// set first row by a vector
//    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));
//
//    /// set second column by a vector
//    M.set_column(1, Vector3D(5.5, 5.5, 5.5));
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
//    /// get the the element at row 1 and column 2
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
//    ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'
//
//    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).


    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;
    //      - recover rotation R and t.

    // to be deleted
    // TODO step 1: Estimate the fundamental matrix F;
    // we apply a translation and a scaling on the image coordinates qi = Tpi, qi'=T'pi'.

    // 1. For every image we calculate the origin of the new coordinate system located at the centroid of the image points.
    Matrix33 T0 = transf_matrix(points_0);
//    std::cout << "T0: " << T0 << std::endl;
    Matrix33 T1 = transf_matrix(points_0);

    // create the structures to hold the NEW NORMALIZED POINTS.
    std::vector<Vector2D> norm_points_0; /// input: 2D image points in the 1st image.
    std::vector<Vector2D> norm_points_1;

    // make the new point qi = T*Pi -- first image
    for (const auto& p0: points_0){
        Vector3D q0_homo = T0 * p0.homogeneous()   ;
        Vector2D q0 = q0_homo.cartesian();
        norm_points_0.emplace_back(q0);
    }
    // make the new point qi = T*Pi -- second image
    for (const auto& p1: points_1){
        Vector3D q1_homo = T1 * p1.homogeneous()   ;
        Vector2D q1 = q1_homo.cartesian();
        norm_points_1.emplace_back(q1);
    }

    /*** ------ Using the Normalized coordinates we will compute the new Fundamentral matrix Fq with the regular 8-point algorithm
    including the rank-2 approximation. However, the matrix Fq is the fundamental matrix for the normalized coordinates.
     We will need to denormalise it, F = ...  ------ ***/

    // Estimate the Fundamental Matrix for Normalized coords Fq, given two images. 8 POINT ALGORITHM
    // QUESTION FOR LIANG : We can only know the Fundamental matrix UP TO A SCALE. ?????????? --> HOW DOES IT AFFECT?
    // Construct W matrix with 24 points
    int m = points_0.size(), n = 9;
    Matrix W(m, n, 0.0);

    for (int i=0; i < m; i++){
        // select the elements
        Vector2D q1 = norm_points_1[i];
        Vector2D q0 = norm_points_0[i];
        W.set_row(i, {q0.x()*q1.x(), q0.x()*q1.x(), q1.x(), q0.x()*q0.y(), q0.y()*q1.y(), q1.y(), q0.x(), q0.y(), 1});
    }
    //std::cout << "W: " << W << std::endl;
    //Compute the SVD decomposition of W
    Matrix U(m, m, 0.0);
    Matrix S(m, n, 0.0); // diagonal matrix that contains the singular values in descending order.
    Matrix V(n, n, 0.0); // stores the corresponding singular vectors. We choose the singular vector p12 which minimises the error.
    svd_decompose(W, U, S, V);

    // PRINTS
    //std::cout << "U: \n" << U << std::endl;
    //std::cout << "S: \n" << S << std::endl;
    //std::cout << "V: \n" << V << std::endl;

    // take the last column of V as the solution for f
    Vector fs = V.get_column(V.cols()-1);
    //std::cout << "fs: \n" << fs << std::endl;

    Matrix33 Fq ;
    Fq.set_row(0, {fs[0], fs[1], fs[2]});
    Fq.set_row(1, {fs[3], fs[4], fs[5]});
    Fq.set_row(2, {fs[6], fs[7], fs[8]});

    //std::cout << "Fq: \n" << Fq << std::endl;

    // Enforce Rank 2 constraint by again using SVD decomposition
    Matrix Uq(3, 3, 0.0);
    Matrix Sq(3, 3, 0.0);
    Matrix Vq(3, 3, 0.0);
    svd_decompose(Fq, Uq, Sq, Vq);

    //std::cout << "Uq: \n" << Uq << std::endl;
    //std::cout << "Sq: \n" << Sq << std::endl;
    //std::cout << "Vq: \n" << Vq << std::endl;

    // enforce rank 2 --> manipulate the diagonal S
    Sq[2][2] = 0;
    // recompute Fq after enforcement
    Fq = Uq * Sq * Vq.transpose();
    //std::cout << "Fq: \n" << Fq << std::endl;


    // Last step: DENORMALIZATION Fq to be F. F = T′TFqT
    Matrix33 F = T1.transpose() * Fq * T0;
    //std::cout << "F: \n" << F << std::endl;

    // TODO: Intermediate step - The recovered F is up to scale. Please scale F such that F(2, 2) = 1.0 after denormalization.
    // so probably take the last element and divide everything with that
    F = F / F[2][2];
//    std::cout << "F: \n" << F << std::endl;


    // TODO step 2: compute the essential matrix E; E = KTFK

    // consruct the K matrix --> same for both cameras
    Matrix33 K(fx, 0, cx, 0, fy, cy, 0, 0, 1);
    //std::cout << "K: \n" << K << std::endl;

    Matrix33 E = K.transpose() * F * K;  // 5 degrees of freedom --> it encodes R, and t (extrinsics).
//    std::cout << "E: \n" << E << std::endl;

    // Decomposition of E into R and t.
    // we define two matrices that we will use in the decomposition of E --> W and Z
    Matrix33 We(0,-1,0,1,0,0,0,0,1);
    Matrix33 Ze(0,1,0,-1,0,0,0,0,0);

    // get the matrices ready for the decomposition of E
    Matrix Ue(3, 3, 0.0);
    Matrix Se(3, 3, 0.0);
    Matrix Ve(3, 3, 0.0);
    svd_decompose(E, Ue, Se, Ve);
//    std::cout << "Ue: \n" << Ue << std::endl;
//    std::cout << "Se: \n" << Se << std::endl;
//    std::cout << "Ve: \n" << Ve << std::endl;


    // to be deleted
    // TODO step 3: recover rotation R and t.
    Matrix33 tr = Ue * Ze * Ue.transpose();
    Vector3D tr1{-tr[1][2], tr[0][2], -tr[0][1]};
    Vector3D tr2 = -tr1;
    // INTERMEDIATE CHECK for t: we can find directly t vector from the U Matrix. t = +- u3 (third column of U)
    Vector3D tr_adele = Ue.get_column(2);

    Matrix33 R1 = Ue * We * Ve.transpose();
    // TODO: check for the determinant of R --> maybe you neglect one this way
    if (!(determinant(R1) < 1.00001 && determinant(R1) > (1 - 0.00001))){
        std::cout << "exit the code --> the determinant of R matrix is wrong. " << std::endl;
        exit(20);
    }
    Matrix33 R2 = Ue * We.transpose() * Ve.transpose();
    if (!(determinant(R2) < 1.00001 && determinant(R2) > (1 - 0.00001))){
        std::cout << "exit the code --> the determinant of R matrix is wrong. " << std::endl;
        exit(20);
    }


    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // create Rt for all possible combinations of R and t.
    Matrix34 R1t1 = construct_Rt(R1, tr1);
    Matrix34 R1t2 = construct_Rt(R1, tr2);
    Matrix34 R2t1 = construct_Rt(R2, tr1);
    Matrix34 R2t2 = construct_Rt(R2, tr2);
    // mutliply K with Rt and find the 3D points
    Matrix34 KR1t1 = K * R1t1;
    Matrix34 KR1t2 = K * R1t2;
    Matrix34 KR2t1 = K * R2t1;
    Matrix34 KR2t2 = K * R2t2;

    //create projection matrix M for the first image
    Matrix33 I = identity(3, 1);
    Vector3D ti{0,0,0};
    Matrix34 M0 = K * construct_Rt(I, ti);
//    std::cout << "M0: " << M0 << std::endl;

    std::vector<Vector3D> points_3d_1 = reconstruct3Dpoints(points_0, points_1, KR1t1,M0);
    std::vector<Vector3D> points_3d_2 = reconstruct3Dpoints(points_0, points_1, KR1t2,M0);
    std::vector<Vector3D> points_3d_3 = reconstruct3Dpoints(points_0, points_1, KR2t1,M0);
    std::vector<Vector3D> points_3d_4 = reconstruct3Dpoints(points_0, points_1, KR2t2,M0);

    // call the function to calculate the number of points in front of both images
    int c1 = count_positive_z(points_3d_1);
    int c2 = count_positive_z(points_3d_2);
    int c3 = count_positive_z(points_3d_3);
    int c4 = count_positive_z(points_3d_4);

    // printing statements to check the number of positive points
//    std::cout << "c1: " << c1 << std::endl;
//    std::cout << "c2: " << c2 << std::endl;
//    std::cout << "c3: " << c3 << std::endl;
//    std::cout << "c4: " << c4 << std::endl;

    // decide on the correct R and t and save it in the variables of the functions.
    if (c1 > c2 && c1 > c3 && c1 > c4){R = R1;t = tr1;points_3d = points_3d_1;}
    if (c2 > c1 && c2 > c3 && c2 > c4){R = R1;t = tr2;points_3d = points_3d_2;}
    if (c3 > c2 && c3 > c1 && c3 > c4){R = R2;t = tr1;points_3d = points_3d_3;}
    if (c4 > c2 && c4 > c3 && c4 > c1){R = R2;t = tr2;points_3d = points_3d_4;}







    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       There are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
//    return points_3d.size() > 0;
    return true;
}

// functions created by the team

// function that returns the centroid of a vector that holds image points
Matrix33 transf_matrix(const std::vector<Vector2D>& im_points){
    unsigned int n0 = im_points.size(); // divide the sum with n and get the centroid
    double p0x = 0, p0y = 0, c0x, c0y;
    for (auto p0: im_points){
        p0x = p0x + p0.x();
        p0y = p0y + p0.y();
    }
    // The origin of our previous CRS is in the left corner so every point will be positive
    // So the translation vector will be negative
    c0x = p0x/n0;
    c0y = p0y/n0;

    // translation vector
    Vector2D t0{-c0x,-c0y};
    // set up the scaling factor --> sqrt(2)/average_dist of the points to the new CRS origin
    double s0, adist0 = 0;
    for (auto p0: im_points){
        adist0 = adist0 + sqrt(pow((p0.x() - c0x), 2) + pow((p0.y() - c0y), 2));
    }
    s0 = sqrt(2)/(adist0/n0);

    // prepare the Transformation matrix (T) for the image
    return {s0, 0, t0.x(), 0, s0, t0.y(), 0, 0, 1};
}

Matrix34 construct_Rt(const Matrix33& R, const Vector3D& t){
    Matrix34 Rt;
    Rt.set_column(0,R.get_column(0));
    Rt.set_column(1,R.get_column(1));
    Rt.set_column(2,R.get_column(2));
    Rt.set_column(3, t);
//    std::cout<< "Rt: " << Rt << std::endl;
    return Rt;
}


std::vector<Vector3D> reconstruct3Dpoints(const std::vector<Vector2D>& img0, const std::vector<Vector2D>& img1, const Matrix34& KRt, const Matrix34& M0){
    // construct Matrix A
    std::vector<Vector3D> points3d;
    const Matrix34& M1=KRt;
    int count = 0;
    for (int i=0; i< img0.size(); i++) {
        // compute the Matrix A
        Matrix44 A;
        A.set_column(0, (img0[i].x() * M0.get_row(2) - M0.get_row(0)));
        A.set_column(0, (img0[i].y() * M0.get_row(2) - M0.get_row(1)));
        A.set_column(0, (img1[i].x() * M1.get_row(2) - M1.get_row(0)));
        A.set_column(0, (img1[i].x() * M1.get_row(2) - M1.get_row(1)));
        Matrix44 U, S, V;
        Matrix Ua(4, 4, 0.0);
        Matrix Sa(4, 4, 0.0);
        Matrix Va(4, 4, 0.0);
        svd_decompose(A, Ua, Sa, Va);
        Vector4D P0 = Va.get_column(3);
        Vector3D p0_3d = P0.cartesian();
        points3d.emplace_back(p0_3d);
    }
    return points3d;
}

int count_positive_z(const std::vector<Vector3D>& points3d){
    int count = 0;
    for (auto p: points3d){
        if (p.z() > 0){
//            // if positive z coord, go and check if the same point in the other image has also a positive coord and add 1 to the counter
//            Vector4D X1 = inverse(KRt) * img1[i].homogeneous();
//            Vector3D p1_3d = X1.cartesian();
//            if (p1_3d.z() > 0){
//                // add 1 to the counter if the corresponding points have positive coords
            count++;
//            }
        }
    }
    return count;
}
