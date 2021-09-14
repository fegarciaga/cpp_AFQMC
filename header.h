#include<iostream>
#include<Eigen/Dense>
#include<unsupported/Eigen/MatrixFunctions>
#include<Eigen/Eigenvalues>
#include<Eigen/StdVector>
#include<cmath>
#include<complex>
#include<random>
#include<vector>

Eigen::MatrixXcf H_K(int Lx, int Ly, int Lz, double kx, double ky, double kz, double tx1, double ty1, double tz1, double tz2);
void initialization(int Lx, int Ly, int Lz, int N_up, int N_dn, int N_wlk, double deltau, double U, double kx, double ky, double kz, double tx1, double ty1, double tz1, double tz2, Eigen::MatrixXcf & H_k, Eigen::MatrixXcf & exp_hk, Eigen::MatrixXcf & Phi_T, std::vector<Eigen::MatrixXcf, Eigen::aligned_allocator<Eigen::MatrixXcf>> & Phi, Eigen::MatrixXd & o, Eigen::MatrixXd & w, Eigen::MatrixXd & aux_field, double fac_norm);
void halfK(Eigen::MatrixXcf & phi, double & w, double & O, Eigen::MatrixXcf & Proj_k_half, Eigen::MatrixXcf Phi_T, Eigen::MatrixXcf & invO_up, Eigen::MatrixXcf & invO_dn, int N_up, int N_dn, int N_sites);
void V(Eigen::MatrixXcf & phi, Eigen::MatrixXcf phi_T, int N_up, int N_dn, double & O, double & w, Eigen::MatrixXcf & invO_matrix_up, Eigen::MatrixXcf & invO_matrix_dn, Eigen::MatrixXd & aux_fld);
double measure(Eigen::MatrixXcf H_k, Eigen::MatrixXcf phi, Eigen::MatrixXcf Phi_T, Eigen::MatrixXcf invO_matrix_up, Eigen::MatrixXcf invO_matrix_dn, int N_up, int N_dn, int N_sites, double U);
void pop_contrl(std::vector<Eigen::MatrixXcf, Eigen::aligned_allocator<Eigen::MatrixXcf>> & Phi, Eigen::MatrixXd & w, Eigen::MatrixXd & o, int N_wlk, int N_sites, int N_par);
void stblz(std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd>> & Phi, int N_wlk, Eigen::MatrixXd & o, int N_up, int N_dn, int N_sites);
void stp_wlk(std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd>> & phi, int N_wlk, int N_sites, Eigen::MatrixXd & w , Eigen::MatrixXd & O, double E, double W, Eigen::MatrixXcf & H_k, Eigen::MatrixXcf & Proj_k_half, bool flag_mea, Eigen::MatrixXd & Phi_T, int N_up, int N_dn, double U, double fac_norm, Eigen::MatrixXd & aux_fld);
void CPMC(int Lx,int Ly,int Lz,int N_up,int N_dn,double kx,double ky,double kz,double U,double tx,double ty,double tz,double ty2,double deltau,int N_wlk, int N_blksteps,int N_eqblk,int N_blk,int itv_modsvd,int itv_pc,int itv_Em);
