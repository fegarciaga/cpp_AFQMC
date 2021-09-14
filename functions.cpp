#include "header.h"

Eigen::MatrixXcf H_K(int Lx, int Ly, int Lz, double kx, double ky, double kz, double tx1, double ty1, double tz1, double tz2)
{
  int r=0;
  int N_sites=Lx*Ly*Lz;
  Eigen::MatrixXcf H_k(N_sites, N_sites);
  H_k.setZero();
  std::complex<double> Kx=(0,1*M_PI*kx);
  std::complex<double> Ky=(0,1*M_PI*ky);
  std::complex<double> Kz=(0,1*M_PI*kz);
  for(int mz=0; mz<Lz;mz++){
    for(int iy=0;iy<Ly;iy++){
      for(int jx=0;jx<Lx;jx++){
        if(Lx!=1){
          if(jx==0){
            H_k(r,r+Lx-1)-=tx1*std::exp(Kx);
            H_k(r,r+1)-=tx1;
          }
          else if(jx==Lx-1){
            H_k(r,r-Lx+1)-=tx1*std::exp(-Kx);
            H_k(r,r-1)-=tx1;
          }
          else{
            H_k(r,r+1)=-tx1;
            H_k(r,r-1)=-tx1;
          }
        }
        if(Ly!=1){
          if(iy==0){
            H_k(r,r+Lx*(Ly-1)-=ty1*std::exp(Ky);
            H_k(r,r+Lx)-=ty1;
          }
          else if(iy==Ly-1){
            H_k(r,r-Lx*(Ly-1))-=ty1*std::exp(-Ky);
            H_k(r,r-Lx)-=ty1;
          }
          else{
            H_k(r,r+Lx)=-ty1;
            H_k(r,r-Lx)=-ty1;
          }
        }
        if(Lz!=1){
          if(mz==0){
            H_k(r,r+Lx*Ly*(Lz-1)-=tz1*std::exp(Kz);
            H_k(r,r+Lx*Ly)-=tz1;
          }
          else if(mz==Lz-1){
            H_k(r,r-Lx*Ly*(Lz-1))-=tz1*std::exp(-Kz);
            H_k(r,r-Lx*Ly)-=tz1;
          }
          else{
            H_k(r,r+Lx*Ly)=-tz1;
            H_k(r,r-Lx*Ly)=-tz1;
          }
        }
        if(Lx>1){
          if(Ly>1){
            if(jx==0){
              if(iy==0){
                H_k(r,r+Lx-1+(Ly-1)*Lx)-=tz2*std::exp(Kx)*std::exp(Ky);
                H_k(r,r+2*Lx-1)-=-tz2*std::exp(Kx);
                H_k(r,r+1+(Ly-1)*Lx)-=tz2*std::exp(Ky);
                H_k(r,r+1+Lx)=-tz2;
              }
              else if(iy==Ly-1){
                H_k(r,r+(Lx-1)-(Ly-1)*Lx)-=tz2*std::exp(Kx)*std::exp(-Ky);
                H_k(r,r+1-(Ly-1)*Lx)-=tz2*std::exp(-Ky);
                H_k(r,r-1)-=tz2*std::exp(Kx);
                H_k(r,r+1-Lx)=-tz2;
              }
              else{
                H_k(r,r+2*Lx-1)-=tz2*std::exp(Kx);
                H_k(r,r-1)-=tz2*std::exp(Kx);
                H_k(r,r+1+Lx)= -tz2;
                H_k(r,r+1-Lx)= -tz2;
              }
            }
            else if (jx==Lx-1){
              if(iy==0){
                H_k(r,r-(Lx-1)+(Ly-1)*Lx)-=tz2*std::exp(-Kx)*exp(Ky);
                H_k(r,r+1)-=tz2*std::exp(-Kx);
                H_k(r,r+(Ly-1)*Lx-1)-=tz2*std::exp(Ky);
                H_k(r,r-1+Lx)=-tz2;
              }
              else if(iy==Ly-1){
                H_k(r,r-(Lx-1)-(Ly-1)*Lx)-=tz2*std::exp(-Kx)*std::exp(-ky);
                H_k(r,r-(2*Lx-1))-= tz2*std::exp(-Kx);
                H_k(r,r-1-(Ly-1)*Lx)-= tz2*std::exp(-Ky);
                H_k(r,r-1-Lx)=-tz2;
              }
              else {
                H_k(r,r+1)-=tz2*std::exp(-Kx);
                H_k(r,r-(2*Lx-1))-=tz2*std::exp(-Kx);
                H_k(r,r-1-Lx)=-tz2;
                H_k(r,r-1+Lx)=-tz2;
              }
            }
            else if(iy==0){
              if((jx+1)%Lx>1){
                H_k(r,r-1+(Ly-1)*Lx)-=tz2*std::exp(Ky);
                H_k(r,r+1+(Ly-1)*Lx)-=tz2*std::exp(Ky);
                H_k(r,r-1+Lx)=-tz2;
                H_k(r,r+1+Lx)=-tz2;
              }
            }
            else if(iy==Ly-1){
              if((jx+1)%Lx>1){
                H_k(r,r-1-(Ly-1)*Lx)-= tz2*std::exp(-Ky);
                H_k(r,r+1-(Ly-1)*Lx)-= tz2*std::exp(-Ky);
                H_k(r,r-1-Lx)=-tz2;
                H_k(r,r+1-Lx)=-tz2;
              }
            }
            else{
              H_k(r,r+1+Lx)=-tz2;
              H_k(r,r-1+Lx)=-tz2;
              H_k(r,r+1-Lx)=-tz2;
              H_k(r,r-1-Lx)=-tz2;
            }
          }
        }
        r+=1;
      }
    }
  }
  return H_k;
}

void initialization(int Lx, int Ly, int Lz, int N_up, int N_dn, int N_wlk, double deltau, double U, double kx, double ky, double kz, double tx1, double ty1, double tz1, double tz2, Eigen::MatrixXcf & H_k, Eigen::MatrixXcf & exp_hk, Eigen::MatrixXcf & Phi_T, std::vector<Eigen::MatrixXcf, Eigen::aligned_allocator<Eigen::MatrixXcf>> & Phi, Eigen::MatrixXd & o, Eigen::MatrixXd & w, Eigen::MatrixXd & aux_field, double fac_norm)
{
  int N_sites=Lx*Ly*Lz;
  int N_par=N_up+n-dn;
  H_k=H_K(Lx, Ly, Lz, kx, ky, kz, tx1, ty1, tz1, tz2);
  exp_hk=(-0.5*deltau*H_k).exp();
  Eigen::EigenSolver<MatrixXcf> es(H_k);
  // Construct trial wave-functions
  Phi_T.block<N_sites,N_up>(0,0)=(es.eigenvectors()).block<N_sites,N_up>(0,0);
  Phi_T.block<N_sites,N_dn>(0,N_up+1)=(es.eigenvectors()).block<N_sites,N_dn>(0,0);
  // Calculate kinetic energy
  double E=0;
  for(int ii=0;ii<N_up;ii++){
      E+=es.eigenvalues()[ii];
  }
  for(int ii=0;ii<N_dn;ii++){
      E+=es.eigenvalues()[ii];
  }
  // Calculate potential energy
  Eigen::MatrixXd n_up;
  Eigen::MatrixXd n_dn;
  n_up=Phi_T.block<N_sites,N_up>(0,0)*(Phi_T.block<N_sites,N_up>(0,0)).adjoint();
  n_dn=Phi_T.block<N_sites,N_dn>(0,N_up)*(Phi_T.block<N_sites,N_dn>(0,N_up)).adjoint();
  for(int ii=0;ii<N_sites,ii++){
    E+=U*n_up(ii,ii)*N_dn(ii,ii);
  }
  // Assemble population of walkers
  for(int ii=0; ii<N_wlk; ii++){
    Phi[ii]=Phi_T;
  }
  w.setOnes();
  o.setOnes();
  // Set auxiliary field
  fac_norm=(std::real(E)-0.5*U*N_par)*deltau;
  double gamma=std::acosh(std::exp(0.5*deltau*U));
  double aux1=-1;
  for(int ii=0;ii<2;ii++){
    aux1*=-1;
    double aux2=-1;
    for(int jj=0;jj<2;jj++){
      aux2*=-1;
      aux_field(ii,jj)=std::exp(gamma*aux1*aux2);
    }
  }
  // Set the random number generator
  std::mt19937 gen(1);
  std::uniform_real_distribution<double> dis(0,1);
  return;
}

void halfK(Eigen::MatrixXcf & phi, double & w, double & O, Eigen::MatrixXcf & Proj_k_half, Eigen::MatrixXcf Phi_T, Eigen::MatrixXcf & invO_up, Eigen::MatrixXcf & invO_dn, int N_up, int N_dn, int N_sites)
{
  phi=Proj_k_half*phi;
  invO_up=((Phi_T.block<N_sites,N_up>(0,0)).adjoint()*phi.block<N_sites,N_up>(0,0)).inverse();
  invO_dn=((Phi_T.block<N_sites,N_dn>(0,N_up)).adjoint()*phi.block<N_sites,n_dn>(0,N_up)).inverse();
  std::complex<double> O_new=1/(invO_up.determinant()*invO_dn.determinant());
  std::complex<double> O_rat=O_new/O;
  if(O_rat<0){
    O=O_new;
    w*=std::real(O_rat);
  }
  else{
    w=0;
  }
  return;
}

void V(Eigen::MatrixXcf & phi, Eigen::MatrixXcf phi_T, int N_up, int N_dn, double & O, double & w, Eigen::MatrixXcf & invO_matrix_up, Eigen::MatrixXcf & invO_matrix_dn, Eigen::MatrixXd & aux_fld)
{
  // Calculate Green's functions
  Eigen::MatrixXcf Gii(2,1);
  Gii.setZero();
  Eigen::MatrixXd RR(2,2);
  RR.setOnes();
  Eigen::MatrixXd matone(2,2);
  matone.setOnes();
  Eigen::MatrixXd G(2,2);
  G.setZero();
  Eigen::MatrixXd temp1_up=phi.block<1,N_up>(0,0)*invO_matrix_up;
  Eigen::MatrixXd temp1_dn=phi.block<1,N_dn>(0,N_up)*invO_matrix_dn;
  Eigen::MatrixXd temp2_up=invO_matrix_up*(phi_T.block<1,N_up>(0,0)).adjoint();
  Eigen::MatrixXd temp2_dn=invO_matrix_dn*(phi_T.block<1,N_dn>(0,N_up)).adjoint();
  Gii(0,0)=temp1_up*(phi_T.block<1,N_up>(0,0)).adjoint();
  Gii(1,0)=temp1_dn*(phi_T.block<1,N_dn>(0,N_up)).adjoint();
  G.col(0)=Gii.col(0);
  G.col(1)=Gii.col(0);
  RR=((aux_fld-matone).cwiseProduct(G))+matone;
  // Implement importance sample
  Eigen::MatrixXd O_ratio=RR.row(0).cwiseProduct(RR.row(1));
  if(std::real(O_ratio(0,0))<0){
    O_ratio(0,0)=0;
  }
  if(std::real(O_ratio(0,1))<0){
    O_ratio(0,1)=0;
  }
  double sum=O_ratio(0,0)+O_ratio(0,1);
  if(sum<=0){
    w=0;
  }
  if(w>0){
    w*=0.5*sum;
    int x_spin=0;
    if(O_raio(0,0)/sum<dis(gen)){
      x_spin=1;
    }
    phi.block<1,N_up>(0,0)=phi.block<1,N_up>(0,0)*aux_fld(0,x_spin);
    phi.block<1,N_dn>(0,N_up)=phi.block<1,N_dn>(0,N_up)*aux_fld(1,x_spin);
    O*=O_ratio(0,x_spin);
    invO_matrix_up=invO_matrix_up+(1-aux_fld(0,x_spin))/RR(0,x_spin)*temp2_up*temp1_up;
    invO_matrix_dn=invO_matrix_dn+(1-aux_fld(1,x_spin))/RR(1,x_spin)*temp2_dn*temp1_dn;
  }
  return;
}

double measure(Eigen::MatrixXcf H_k, Eigen::MatrixXcf phi, Eigen::MatrixXcf Phi_T, Eigen::MatrixXcf invO_matrix_up, Eigen::MatrixXcf invO_matrix_dn, int N_up, int N_dn, int N_sites, double U)
{
  double E=0;
  Eigen::MatrixXcf G_up=(phi.block<N_sites,N_up>(0,0)*invO_matrix_up)*(phi_T.block<N_sites,N_up>(0,0)).adjoint();
  Eigen::MatrixXcf G_dn=(phi.block<N_sites,N_dn>(0,N_up)*invO_matrix_up)*(phi_T.block<N_sites,N_dn>(0,N_up)).adjoint();
  E+=U*((G_up.diagonal()).transpose()*G_dn.diagonal());
  E+=(H_k.cwiseProduct(G_up+G_dn)).sum();
  return E;
}

void pop_contrl(std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd>> & Phi, Eigen::MatrixXd & w, Eigen::MatrixXd & o, int N_wlk, int N_sites, int N_par)
{
  std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd>>  new_Phi;
  Eigen::MatrixXd new_O(N_wlk,1);
  double d=N_wlk/w.sum();
  double sum_w=-dis(gen);
  int n_wlk=0;
  for(int ii=0; ii<N_wlk; ii++){
    sum_w+=w(ii,0)*d;
    int n=std::ceil(sum_w);
    for(int jj=n_wlk; jj<n;jj++){
      new_Phi[jj]=Phi[ii];
      new_O(jj,0)=o(ii,0);
    }
  }
  Phi=new_Phi;
  o=new_O;
  w.setOnes();
  return;
}

void stblz(std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd>> & Phi, int N_wlk, Eigen::MatrixXd & o, int N_up, int N_dn, int N_sites)
{
  for(int ii=0; ii<N_wlk;ii++){
    Eigen::MatrixXcf thinQ((Phi[ii].block<N_sites,N_up>(0,0)).rows(),(Phi[ii].block<N_sites,N_up>(0,0)).cols()), q((Phi[ii].block<N_sites,N_up>(0,0)).rows(),(Phi[ii].block<N_sites,N_up>(0,0)).cols());
    Eigen::HouseholderQR<Eigen::MatrixXcf> householderQR(Phi[ii].block<N_sites,N_up>(0,0));
    q=householderQR.householderQ();
    thinQ.setIdentity();
    Q=householderQR.householderQ()*thinQ;
    R=Q.transpose()*Phi[ii].block<N_sites,N_up>(0,0);
    Phi[ii].block<N_sites,N_up>(0,0)=Q;
    o(ii,0)/=R.determinant();
    Eigen::MatrixXcf thinQ((Phi[ii].block<N_sites,N_dn>(0,N_up)).rows(),(Phi[ii].block<N_sites,N_dn>(0,N_up)).cols()), q((Phi[ii].block<N_sites,N_dn>(0,N_up)).rows(),(Phi[ii].block<N_sites,N_dn>(0,N_up)).cols());
    Eigen::HouseholderQR<Eigen::MatrixXcf> householderQR(Phi[ii].block<N_sites,N_dn>(0,N_up));
    q=householderQR.householderQ();
    thinQ.setIdentity();
    Q=householderQR.householderQ()*thinQ;
    R=Q.transpose()*Phi[ii].block<N_sites,N_dn>(0,N_up);
    Phi[ii].block<N_sites,N_dn>(0,N_up)=Q;
    o(ii,0)/=R.determinant();
  }
  return;
}

void stp_wlk(std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd>> & phi, int N_wlk, int N_sites, Eigen::MatrixXd & w , Eigen::MatrixXd & O, double E, double W, Eigen::MatrixXcf & H_k, Eigen::MatrixXcf & Proj_k_half, bool flag_mea, Eigen::MatrixXd & Phi_T, int N_up, int N_dn, double U, double fac_norm, Eigen::MatrixXd & aux_fld)
{
  Eigen::MatrixXd e(N_wlk,1);
  e=setZero();
  Eigen::MatrixXcf invO_up(N_up,N_up);
  invO_up.setZero();
  Eigen::MatrixXcf invO_dn(N_dn,N_dn);
  invO_dn.setZero();
  for(int ii=0; ii<N_wlk; ii++){
    if(w(ii,0)>0){
      w(ii,0)*=std::exp(fac_norm);
      halfK(phi[ii], w(ii,0), O(ii,0), Proj_k_half, Phi_T, invO_up, invO_dn, N_up, N_dn, N_sites);
      if(w(ii,0)>0){
        for(int jj=0;jj<N_sites;jj++){
          if(w(ii,0)>0){
            V(phi[ii].row(jj), phi_T.row(jj), N_up, N_dn, O(ii,0), w(ii,0), invO_up, invO_dn, aux_fld);
          }
        }
      }
      if(w(ii,0)>0){
        halfK(phi[ii], w(ii,0), O(ii,0), Proj_k_half, Phi_T, invO_up, invO_dn, N_up, N_dn, N_sites);
        if(w(ii,0)>0){
          if (flag_mea==true){
            e(ii,0)=measure(H_k, phi[ii], Phi_T, invO_up, invO_dn, N_up, N_dn, N_sites, U);
          }
        }
      }
  }
  if(flag_mea==true){
    for(int ii=0; ii<N_wlk; ii++){
      if(w(ii,0)>0){
        E+=e(ii,0)*w(ii,0);
        W+=w(ii,0);
      }
    }
  }
}

void CPMC(int Lx,int Ly,int Lz,int N_up,int N_dn,double kx,double ky,double kz,double U,double tx,double ty,double tz,double ty2,double deltau,int N_wlk, int N_blksteps,int N_eqblk,int N_blk,int itv_modsvd,int itv_pc,int itv_Em)
{
  int N_par=N_up+N_dn;
  double fac_norm=0;
  Eigen::MatrixXd aux_field(2,2);
  Eigen::MatrixXcf H_k(N_sites, N_sites);
  Eigen::MatrixXcf Phi_T(N_sites, N_par);
  Eigen::MatrixXd w(N_wlk,1);
  Eigen::MatrixXd o(N_wlk,1);
  std::vector<Eigen::MatrixXcf, Eigen::aligned_allocator<Eigen::MatrixXcf>> Phi(N_wlk);
  bool=flag_mea=false;
  initialization(Lx, Ly, Lz, N_up, N_dn, N_wlk,deltau, U, kx, ky, kz, tx1, ty1, tz1, tz2, H_k, exp_hk, Phi_T, Phi, o, w, aux_field, fac_norm);
  double E=0;
  double W=0;
  Eigen::MatrixXcf E_blk(N_blk,1);
  Eigen::MatrixXcf W_blk(N_blk,1);
  for(int ii=0; ii<N_eqblk;ii++){
    for(int jj=0; jj<N_blksteps;jj++){
      stp_wlk(Phi, N_wlk, N_sites, w , o, E, W, H_k, Proj_k_half, flag_mea, Phi_T, N_up, N_dn, U, fac_norm, aux_fld);
      if((jj+1)%itv_modsvd==0){
        stblz(Phi, N_wlk, o, N_up, N_dn, N_sites);
      }
      if((jj+1)%itv_pc==0){
        pop_contrl(Phi, w, o, N_wlk, N_sites, N_par);
      }
    }
  }
  for(int ii=0; ii<N_blk;ii++){
    for(int jj=0; jj<N_blksteps;jj++){
      if ((jj+1)%itv_Em==0){
        flag_mea=true;
      }
      else{
        flag_mea=false;
      }
      stp_wlk(Phi, N_wlk, N_sites, w , o, E_blk(ii,0), W_blk(ii,0), H_k, Proj_k_half, flag_mea, Phi_T, N_up, N_dn, U, fac_norm, aux_fld);
      if((jj+1)%itv_modsvd==0){
        stblz(Phi, N_wlk, o, N_up, N_dn, N_sites);
      }
      if((jj+1)%itv_pc==0){
        pop_contrl(Phi, w, o, N_wlk, N_sites, N_par);
      }
      if((jj+1)%itv_Em==0){
        fac_norm=(std::real(E_blk(i_blk,0)/W_blk(i_blk,0))-0.5*U*N_par)*deltau;
      }
    }
    E_blk(i_blk,0)/=W_blk(i_blk,0);
  }
  std::cout<<std::real(E_blk.mean())<<"\t"<<std::sqrt(std::real(E_blk-E_blk.mean()).square().sum()/(E_blk.size()-1))/std::sqrt(N_wlk)<<"\n";
}
