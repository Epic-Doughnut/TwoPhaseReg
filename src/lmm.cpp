#include <RcppEigen.h>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include <vector>
#include "utility.h"

using namespace Rcpp;

typedef Map<VectorXd> MapVecd;
typedef Map<VectorXi> MapVeci;
typedef Map<MatrixXd> MapMatd;
typedef Map<MatrixXi> MapMati;

double LMM_GeneralSplineProfile (MatrixXd& pB, RowVectorXd& p_col_sum, VectorXd& q_row_sum, VectorXd& vcparams, 
	VectorXd& vcparams0, MatrixXd& p, MatrixXd& p0, MatrixXd& P_theta, MatrixXd& q, MatrixXd& logp, 
	VectorXd& Z_theta, VectorXd& X_uni_theta, VectorXd& TX_uni_theta, VectorXd& resi2, MatrixXd& resi, 
	MatrixXd& D, vector<MatrixXd>& V, vector<MatrixXd>& V_inv, vector<VectorXd>& V_inv_1, vector<VectorXd>& V_inv_T, 
	VectorXd& V_inv_11, VectorXd& V_inv_1T, VectorXd& V_inv_TT, VectorXd& L_var, MatrixXd& I_var, 
	VectorXd& TZ_theta, 
	const VectorXd& theta, const VectorXd& vc_initial, const VectorXd& Y, const VectorXd& T, 
	const MatrixXd& U, const MatrixXi& index_obs, const MatrixXd& Bspline_uni, 
	const MatrixXd& Z, const MatrixXd& X_uni, const VectorXi& X_uni_ind, 
	const VectorXi& Bspline_uni_ind, const MatrixXd& p_static, const int& n, const int& n2, const int& m, 
	const int& s, const int& n_minus_n2, const int& nobs, const int& nobs2, const int& X_nc, const int& Z_nc, const int& nT, const int& MAX_ITER, 
	const double& TOL, const int& nZ, const int& nXT, const bool& ZT_flag, const int& ZT_nc, const MatrixXd& ZT)
{
/**** temporary variables **********************************************************************************************************************/
	double tol, V_inv_1R, V_inv_TR, V_inv_R2;
	int iter, idx, idxx;
// /* RT's test code */
// time_t t1, // t2;
// /* RT's test code end */
/**** temporary variables **********************************************************************************************************************/

/**** EM algorithm *****************************************************************************************************************************/

/**** parameter initialization *****************************************************************************************************************/
	vcparams = vc_initial;
	vcparams0 = vc_initial;

	p_col_sum = p_static.colwise().sum();
	for (int j=0; j<s; ++j)
	{
		p.col(j) = p_static.col(j)/p_col_sum(j);
	}
	p0 = p;
	MatrixXd Bspline_Mat(n_minus_n2,m);
 	MatrixXd pthetaOverQ(P_theta.rows(), P_theta.cols());

/**** parameter initialization *****************************************************************************************************************/

/**** calculate residual ***********************************************************************************************************************/
	X_uni_theta = X_uni*theta.segment(1, X_nc);
	Z_theta = Z*theta.segment(nZ, Z_nc);
	TX_uni_theta = X_uni*theta.segment(nXT, X_nc);
	if (ZT_flag)
	{
		TZ_theta = ZT*theta.tail(ZT_nc);
	}
	for (int i=0; i<n2; ++i)
	{
		resi2.segment(index_obs(i,0), index_obs(i, 1)) = Y.segment(index_obs(i,0), index_obs(i, 1));
		resi2.segment(index_obs(i,0), index_obs(i, 1)).array() -= theta(0)+X_uni_theta(X_uni_ind(i))+Z_theta(i);
		resi2.segment(index_obs(i,0), index_obs(i, 1)).noalias() -= T.segment(index_obs(i,0), index_obs(i, 1))*theta(nT);
		resi2.segment(index_obs(i,0), index_obs(i, 1)).noalias() -= T.segment(index_obs(i,0), index_obs(i, 1))*TX_uni_theta(X_uni_ind(i));
		if (ZT_flag)
		{
			resi2.segment(index_obs(i,0), index_obs(i, 1)).noalias() -= T.segment(index_obs(i,0), index_obs(i, 1))*TZ_theta(i);
		}
	}
	for (int i=0; i<n_minus_n2; ++i)
	{
		idx = i+n2;
		idxx = index_obs(idx,0)-nobs2;
		resi.block(idxx, 0, index_obs(idx, 1), 1) = Y.segment(index_obs(idx, 0), index_obs(idx, 1));
		resi.block(idxx, 0, index_obs(idx, 1), 1).array() -= theta(0)+Z_theta(idx);
		resi.block(idxx, 0, index_obs(idx, 1), 1).noalias() -= T.segment(index_obs(idx, 0), index_obs(idx, 1))*theta(nT);
		if (ZT_flag)
		{
			resi.block(idxx, 0, index_obs(idx, 1), 1).noalias() -= T.segment(index_obs(idx, 0), index_obs(idx, 1))*TZ_theta(idx);
		}
		for (int k=1; k<m; ++k)
		{
			resi.block(idxx, k, index_obs(idx, 1), 1) = resi.block(idxx, 0, index_obs(idx, 1), 1);
		}
		for (int k=0; k<m; ++k)
		{
			resi.block(idxx, k, index_obs(idx, 1), 1).array() -= X_uni_theta(k);
			resi.block(idxx, k, index_obs(idx, 1), 1).noalias() -= T.segment(index_obs(idx, 0), index_obs(idx, 1))*TX_uni_theta(k);
		}
	}
/**** calculate residual ***********************************************************************************************************************/

	for (iter=0; iter<MAX_ITER; ++iter)
	{
		// /* RT's test code */
		// time(&t1);
		// /* RT's test code end */

		/**** E-step *******************************************************************************************************************************/

		/**** update pB ****************************************************************************************************************************/
		pB = Bspline_uni*p.transpose();
		/**** update pB ****************************************************************************************************************************/

		/**** update D *****************************************************************************************************************************/
		D(0,0) = vcparams(0);
		D(0,1) = D(1,0) = vcparams(1);
		D(1,1) = vcparams(2);
		/**** update D *****************************************************************************************************************************/

		/**** update V *****************************************************************************************************************************/
		for (int i=0; i<n; ++i)
		{
			V[i].setIdentity();
			V[i] *= vcparams(3);
			V[i].noalias() += U.block(index_obs(i,0), 0, index_obs(i, 1), 2)*D*U.block(index_obs(i, 0), 0, index_obs(i, 1), 2).transpose();
			V_inv[i] = V[i].selfadjointView<Eigen::Upper>().ldlt().solve(MatrixXd::Identity(index_obs(i,1), index_obs(i, 1)));
		}
		/**** update V *****************************************************************************************************************************/

		/**** update P_theta ***********************************************************************************************************************/
		for (int i=0; i<n_minus_n2; ++i)
		{
			idx = i+n2;
			idxx = index_obs(idx,0)-nobs2;
			for (int k=0; k<m; ++k)
			{
				P_theta(i,k) = (resi.block(idxx, k, index_obs(idx, 1), 1).transpose()*V_inv[idx]*resi.block(idxx, k, index_obs(idx, 1), 1))(0, 0);
			}
		}
		P_theta /= -2.;
		P_theta = P_theta.array().exp();
/**** update P_theta ***********************************************************************************************************************/

/**** update q, q_row_sum ******************************************************************************************************************/
		// for (int i=0; i<n_minus_n2; ++i)
		// {
		// 	for (int k=0; k<m; ++k)
		// 	{
		// 		q(i,k) = P_theta(i,k)*pB(Bspline_uni_ind(i+n2),k);
		// 	}
		// }
		for (int i = 0; i < n_minus_n2; ++i)
		{
			Bspline_Mat.row(i) = pB.row(Bspline_uni_ind(i + n2));
		}
		q = P_theta.array() * Bspline_Mat.array();

		q_row_sum = q.rowwise().sum();
		for (int i=0; i<n_minus_n2; ++i)
		{
			q.row(i) /= q_row_sum(i);
		}
/**** update q, q_row_sum ******************************************************************************************************************/

/**** E-step *******************************************************************************************************************************/


/**** M-step *******************************************************************************************************************************/

/**** update variance components parameters ************************************************************************************************/
		for (int i=0; i<n; ++i)
		{
			V_inv_1[i]  = V_inv[i].rowwise().sum();
			V_inv_T[i]  = V_inv[i]*T.segment(index_obs(i,0), index_obs(i, 1));
			V_inv_11(i) = V_inv_1[i].sum();
			V_inv_1T(i) = V_inv_T[i].sum();
			V_inv_TT(i) = V_inv_T[i].dot(T.segment(index_obs(i,0), index_obs(i, 1)));
		}

		L_var.setZero();
		I_var.setZero();

		for (int i=0; i<n2; ++i)
		{
			V_inv_1R 	 = V_inv_1[i].dot(resi2.segment(index_obs(i,0), index_obs(i, 1)));
			V_inv_TR 	 = V_inv_T[i].dot(resi2.segment(index_obs(i,0), index_obs(i, 1)));
			V_inv_R2 	 = (V_inv[i]*resi2.segment(index_obs(i,0), index_obs(i, 1))).squaredNorm();

			L_var(0) 	+= V_inv_1R*V_inv_1R-V_inv_11(i);
			L_var(1) 	+= 2.*(V_inv_TR*V_inv_1R-V_inv_1T(i));
			L_var(2) 	+= V_inv_TR*V_inv_TR-V_inv_TT(i);
			L_var(3) 	+= V_inv_R2-V_inv[i].trace();

			I_var(0,0) 	+= V_inv_11(i)*V_inv_11(i);
			I_var(0,1) 	+= 2.*V_inv_1T(i)*V_inv_11(i);
			I_var(0,2) 	+= V_inv_1T(i)*V_inv_1T(i);
			I_var(0,3) 	+= V_inv_1[i].squaredNorm();
			I_var(1,1) 	+= 2*(V_inv_1T(i)*V_inv_1T(i)+V_inv_11(i)*V_inv_TT(i));
			I_var(1,2) 	+= 2*V_inv_1T(i)*V_inv_TT(i);
			I_var(1,3) 	+= 2*V_inv_1[i].dot(V_inv_T[i]);
			I_var(2,2) 	+= V_inv_TT(i)*V_inv_TT(i);
			I_var(2,3) 	+= V_inv_T[i].squaredNorm();
			I_var(3,3) 	+= (V_inv[i]*V_inv[i]).trace();
		}

		for (int i=0; i<n_minus_n2; ++i)
		{
			idx = i+n2;
			idxx = index_obs(idx,0)-nobs2;

			for (int k=0; k<m; ++k)
			{
				V_inv_1R = (V_inv_1[idx].transpose()*resi.block(idxx, k, index_obs(idx, 1), 1))(0, 0);
				V_inv_TR = (V_inv_T[idx].transpose()*resi.block(idxx, k, index_obs(idx, 1), 1))(0, 0);
				V_inv_R2 = (V_inv[idx]*resi.block(idxx, k, index_obs(idx, 1), 1)).squaredNorm();

				L_var(0) += q(i,k)*V_inv_1R*V_inv_1R;
				L_var(1) += 2.*q(i,k)*V_inv_TR*V_inv_1R;
				L_var(2) += q(i,k)*V_inv_TR*V_inv_TR;
				L_var(3) += q(i,k)*V_inv_R2;
			}

			L_var(0) -= V_inv_11(idx);
			L_var(1) -= 2.*V_inv_1T(idx);
			L_var(2) -= V_inv_TT(idx);
			L_var(3) -= V_inv[idx].trace();

			I_var(0,0) 	+= V_inv_11(idx)*V_inv_11(idx);
			I_var(0,1) 	+= 2.*V_inv_1T(idx)*V_inv_11(idx);
			I_var(0,2) 	+= V_inv_1T(idx)*V_inv_1T(idx);
			I_var(0,3) 	+= V_inv_1[idx].squaredNorm();
			I_var(1,1) 	+= 2*(V_inv_1T(idx)*V_inv_1T(idx)+V_inv_11(idx)*V_inv_TT(idx));
			I_var(1,2) 	+= 2*V_inv_1T(idx)*V_inv_TT(idx);
			I_var(1,3) 	+= 2*V_inv_1[idx].dot(V_inv_T[idx]);
			I_var(2,2) 	+= V_inv_TT(idx)*V_inv_TT(idx);
			I_var(2,3) 	+= V_inv_T[idx].squaredNorm();
			I_var(3,3) 	+= (V_inv[idx]*V_inv[idx]).trace();
		}

		vcparams = I_var.selfadjointView<Eigen::Upper>().ldlt().solve(L_var);
		vcparams += vcparams0;
/**** update variance components parameters ************************************************************************************************/

/**** update p *****************************************************************************************************************************/
		p.setZero();
		pthetaOverQ = P_theta.array().colwise() / q_row_sum.array();

		for (int i=0; i<n_minus_n2; ++i)
		{
			// for (int k=0; k<m; ++k)
			// {
			// 	for (int j=0; j<s; ++j)
			// 	{
			// 		p(k,j) 	+= Bspline_uni(Bspline_uni_ind(i+n2),j)*P_theta(i,k)/q_row_sum(i);
			// 	}
			// }
			p += pthetaOverQ.row(i).transpose() * Bspline_uni.row(Bspline_uni_ind(i + n2));

		}
		p = p.array()*p0.array();
		p += p_static;
		p_col_sum = p.colwise().sum();
		for (int j=0; j<s; ++j)
		{
			p.col(j) /= p_col_sum(j);
		}
/**** update p *****************************************************************************************************************************/

/**** M-step *******************************************************************************************************************************/


/**** calculate the sum of absolute differences between estimates in the current and previous iterations ***********************************/
		tol = (vcparams-vcparams0).array().abs().sum();
		tol += (p-p0).array().abs().sum();
/**** calculate the sum of absolute differences between estimates in the current and previous iterations ***********************************/

/**** update parameters ********************************************************************************************************************/
		vcparams0 = vcparams;
		p0 = p;
/**** update parameters ********************************************************************************************************************/

/**** check convergence ********************************************************************************************************************/
		if (tol < TOL)
		{
			break;
		}
/**** check convergence ********************************************************************************************************************/

// /* RT's test code */
// time(&t2);
// Rcout << iter << '\t' << difftime(t2, t1) << '\t' << tol << endl;
// /* RT's test code end */
	}
/**** EM algorithm *****************************************************************************************************************************/

	if (iter == MAX_ITER)
	{
		return -999.;
	}
	else
	{
/**** calculate the likelihood *************************************************************************************************************/
		double loglik;

		logp = p.array().log();
		for (int k=0; k<m; ++k)
		{
			for (int j=0; j<s; ++j)
			{
				if (p(k,j) <= 0.)
				{
					logp(k,j) = 0.;
				}
			}
		}
		pB = Bspline_uni*logp.transpose();

		loglik = 0.;
		for (int i=0; i<n2; ++i) 
		{
			loglik += pB(Bspline_uni_ind(i),X_uni_ind(i));
		}

		pB = Bspline_uni*p.transpose();

		D(0,0) = vcparams(0);
		D(0,1) = D(1,0) = vcparams(1);
		D(1,1) = vcparams(2);

		for (int i=0; i<n; ++i)
		{
			V[i].setIdentity();
			V[i] *= vcparams(3);
			V[i].noalias() += U.block(index_obs(i,0), 0, index_obs(i, 1), 2)*D*U.block(index_obs(i, 0), 0, index_obs(i, 1), 2).transpose();
			V_inv[i] = V[i].selfadjointView<Eigen::Upper>().ldlt().solve(MatrixXd::Identity(index_obs(i,1), index_obs(i, 1)));
		}

		for (int i=0; i<n_minus_n2; ++i)
		{
			idx = i+n2;
			idxx = index_obs(idx,0)-nobs2;
			for (int k=0; k<m; ++k)
			{
				P_theta(i,k) = (resi.block(idxx, k, index_obs(idx, 1), 1).transpose()*V_inv[idx]*resi.block(idxx, k, index_obs(idx, 1), 1))(0, 0);
			}
		}
		P_theta /= -2.;
		P_theta = P_theta.array().exp();

		// for (int i=0; i<n_minus_n2; ++i)
		// {
		// 	for (int k=0; k<m; ++k)
		// 	{
		// 		q(i,k) = P_theta(i,k)*pB(Bspline_uni_ind(i+n2),k);
		// 	}
		// }
		for (int i = 0; i < n_minus_n2; ++i)
		{
			Bspline_Mat.row(i) = pB.row(Bspline_uni_ind(i+n2));
		}
		q = P_theta.array() * Bspline_Mat.array();
		q_row_sum = q.rowwise().sum();

		loglik += q_row_sum.array().log().sum();
		loglik -= log(2.*M_PI)*nobs/2.;
		for (int i=0; i<n; ++i)
		{
			loglik -= log(V[i].determinant())/2.;
		}
		for (int i=0; i<n2; ++i)
		{
			loglik -= (resi2.segment(index_obs(i,0), index_obs(i, 1)).transpose()*V_inv[i]*resi2.segment(index_obs(i,0), index_obs(i, 1)))(0,0)/2.;
		}
/**** calculate the likelihood *************************************************************************************************************/

		return loglik;
	}
} // LMM_GeneralSplineProfile

// [[Rcpp::export]]
List LMM_GeneralSpline (const MapVecd& Y,
	const MapVecd& T,
	const MapMatd& X,
	const MapMatd& Z,
	const MapVeci& ZT_colid,
	const MapMatd& Bspline,
	const MapMati& index_obs,
	const MapVecd& coef_initial,
	const MapVecd& vc_initial,
	const double& hn,
	const int& MAX_ITER,
	const double& TOL,
	const int& noSE)
{
	/*#############################################################################################################################################*/
	/**** pass arguments from R to cpp *************************************************************************************************************/
	// const MapVecd Y(as<MapVecd>(Y_R));
	// const MapVecd T(as<MapVecd>(T_R));
	// const MapMatd X(as<MapMatd>(X_R));
	// const MapMatd Z(as<MapMatd>(Z_R));
	// const MapVeci ZT_colid(as<MapVeci>(ZT_colid_R));
	// const MapMatd Bspline(as<MapMatd>(Bspline_R));
	// const MapMati index_obs(as<MapMati>(index_obs_R));
	// const MapVecd coef_initial(as<MapVecd>(coef_initial_R));
	// const MapVecd vc_initial(as<MapVecd>(vc_initial_R));
	// const double hn = NumericVector(hn_R)[0];
	// const int MAX_ITER = IntegerVector(MAX_ITER_R)[0];
	// const double TOL = NumericVector(TOL_R)[0];
	// const int noSE = IntegerVector(noSE_R)[0];
	/**** pass arguments from R to cpp *************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** some useful constants ********************************************************************************************************************/
	const int n = Z.rows();  // number of subjects in the first phase
	const int n2 = X.rows(); // number of subjects in the second phase
	const int n_minus_n2 = n-n2; // number of subjects not selected in the second phase
	const int nobs = Y.size(); // number of observations in the first phase
	const int nobs2 = index_obs(n2,0); // number of observations in the second phase
	const int nobs_minus_nobs2 = nobs-nobs2; // number of observations not selected in the second phase
	const int Z_nc = Z.cols(); // number of inexpensive covariates
	const int X_nc = X.cols(); // number of expensive covariates X
	const int nZ = 1+X_nc;
	const int nT = nZ+Z_nc;
	const int nXT = nT+1;
	const int nZT = nXT+X_nc;
	int ncov = nZT;
	bool ZT_flag = false;
	int ZT_nc = 0;
	MatrixXd ZT;
	if (ZT_colid(0) != -1)
	{
		ZT_flag = true;
		ZT_nc = ZT_colid.size();
		ncov += ZT_nc;
		ZT.resize(n,ZT_nc);
		for (int k=0; k<ZT_nc; ++k) {
			ZT.col(k) = Z.col(ZT_colid(k));
		}
	}
	const int s = Bspline.cols(); // number of B-spline functions
	/**** some useful constants ********************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** summarize observed distinct rows of X ****************************************************************************************************/
	// m: number of distinct rows of X
	// X_uni_ind(n2): row index of rows of X in X_uni
	// X_uni(m,X_nc): distinct rows of X
	// X_uni_n(m): count of appearances of each distinct row of X
	VectorXi X_index =indexx_Matrix_Row(X);
	int m = Num_Uni_Matrix_Row(X,X_index);

	VectorXi X_uni_ind(n2);
	MatrixXd X_uni(m,X_nc);
	VectorXi X_uni_n(m);
	Create_Uni_Matrix_Row(X,X_index,X_uni,X_uni_ind,X_uni_n);
	/**** summarize observed distinct rows of X ****************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** summarize observed distinct rows of Bspline **********************************************************************************************/
	// m_B: number of distinct rows of Bspline
	// Bspline_uni_ind(n): row index of rows of Bspline in Bspline_uni
	// Bspline_uni(m_B,s): distinct rows of Bspline
	// Bspline_uni_n(m_B): count of appearances of each distinct row of Bspline
	VectorXi Bspline_index = indexx_Matrix_Row(Bspline);
	int m_B = Num_Uni_Matrix_Row(Bspline,Bspline_index);
	VectorXi Bspline_uni_ind(n);
	MatrixXd Bspline_uni(m_B,s);
	VectorXi Bspline_uni_n(m_B);
	Create_Uni_Matrix_Row(Bspline,Bspline_index,Bspline_uni,Bspline_uni_ind,Bspline_uni_n);
	/**** summarize observed distinct rows of Bspline **********************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** some fixed quantities in the EM algorithm ************************************************************************************************/
	// p
	MatrixXd p_static(m,s);
	p_static.setZero();
	for (int i=0; i<n2; ++i)
	{
		p_static.row(X_uni_ind(i)) 	+= Bspline.row(i);
	}

	// U: random effects design matrix
	MatrixXd U(nobs,2);
	U.col(0).setOnes();
	U.col(1) = T;
	/**** some fixed quantities in the EM algorithm ************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** output ***********************************************************************************************************************************/
	VectorXd theta(ncov); // regression coefficients
	MatrixXd cov_theta(ncov, ncov); // covariance matrix
	bool flag_nonconvergence; // flag of none convergence in the estimation of regression coefficients
	bool flag_nonconvergence_cov; // flag of none convergence in the estimation of covariance matrix
	VectorXd vcparams(4); // variance components parameters
	/**** output ***********************************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** temporary variables **********************************************************************************************************************/
	VectorXd XinvVY(ncov);
	MatrixXd XinvVX(ncov,ncov);
	VectorXd theta0(ncov);
	VectorXd vcparams0(4);
	MatrixXd p(m,s);
	MatrixXd p0(m,s);
	RowVectorXd p_col_sum(s);
	MatrixXd q(n_minus_n2,m);
	VectorXd q_row_sum(n_minus_n2);
	MatrixXd pB(m_B,m);
	MatrixXd P_theta(n_minus_n2,m);
	MatrixXd resi(nobs_minus_nobs2,m);
	VectorXd resi2(nobs2);
	vector<MatrixXd> V(n);
	vector<MatrixXd> V_inv(n);
	vector<VectorXd> V_inv_1(n);
	vector<VectorXd> V_inv_T(n);
	for (int i=0; i<n; ++i) 
	{
		V[i].resize(index_obs(i,1),index_obs(i,1));
		V_inv[i].resize(index_obs(i,1),index_obs(i,1));
		V_inv_1.resize(index_obs(i,1));
		V_inv_T.resize(index_obs(i,1));
	}
	VectorXd V_inv_11(n);
	VectorXd V_inv_1T(n);
	VectorXd V_inv_TT(n);
	MatrixXd D(2,2);
	VectorXd Z_theta(n);
	VectorXd X_uni_theta(m);
	VectorXd TX_uni_theta(m);
	VectorXd TZ_theta;
	if (ZT_flag)
	{
		TZ_theta.resize(n);
	}
	VectorXd L_var(4);
	MatrixXd I_var(4,4);
	double tol,V_inv_1Y,V_inv_TY,V_inv_1R,V_inv_TR,V_inv_R2;
	int iter,idx,idxx;
	// /* RT's test code */
	// time_t t1,t2;
	// /* RT's test code end */
	/**** temporary variables **********************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** EM algorithm *****************************************************************************************************************************/

	/**** parameter initialization *****************************************************************************************************************/
	theta = coef_initial;
	theta0 = coef_initial;
	vcparams = vc_initial;
	vcparams0 = vc_initial;

	p_col_sum = p_static.colwise().sum();
	for (int j=0; j<s; ++j)
	{
		p.col(j) = p_static.col(j)/p_col_sum(j);
	}
	p0 = p;

	flag_nonconvergence = false;
	flag_nonconvergence_cov = false;
 	MatrixXd pthetaOverQ(P_theta.rows(), P_theta.cols());
 	MatrixXd Bspline_Mat(n_minus_n2,m);
	/**** parameter initialization *****************************************************************************************************************/

	/**** residual initialization ******************************************************************************************************************/
	X_uni_theta = X_uni*theta.segment(1,X_nc);
	Z_theta = Z*theta.segment(nZ,Z_nc);
	TX_uni_theta = X_uni*theta.segment(nXT,X_nc);
	if (ZT_flag)
	{
		TZ_theta = ZT*theta.tail(ZT_nc);
	}
	for (int i=0; i<n_minus_n2; ++i)
	{
		idx = i+n2;
		idxx = index_obs(idx,0)-nobs2;
		resi.block(idxx,0,index_obs(idx,1),1) = Y.segment(index_obs(idx,0),index_obs(idx,1));
		resi.block(idxx,0,index_obs(idx,1),1).array() -= theta(0)+Z_theta(idx);
		resi.block(idxx,0,index_obs(idx,1),1).noalias() -= T.segment(index_obs(idx,0),index_obs(idx,1))*theta(nT);
		if (ZT_flag)
		{
			resi.block(idxx,0,index_obs(idx,1),1).noalias() -= T.segment(index_obs(idx,0),index_obs(idx,1))*TZ_theta(idx);
		}
		for (int k=1; k<m; ++k)
		{
			resi.block(idxx,k,index_obs(idx,1),1) = resi.block(idxx,0,index_obs(idx,1),1);
		}
		for (int k=0; k<m; ++k)
		{
			resi.block(idxx,k,index_obs(idx,1),1).array() -= X_uni_theta(k);
			resi.block(idxx,k,index_obs(idx,1),1).noalias() -= T.segment(index_obs(idx,0),index_obs(idx,1))*TX_uni_theta(k);
		}
	}
	/**** residual initialization ******************************************************************************************************************/

	for (iter=0; iter<MAX_ITER; ++iter)
	{
	// /* RT's test code */
	// time(&t1);
	// /* RT's test code end */

	/**** E-step *******************************************************************************************************************************/

	/**** update pB ****************************************************************************************************************************/
		pB = Bspline_uni*p.transpose();
	/**** update pB ****************************************************************************************************************************/

	/**** update D *****************************************************************************************************************************/
		D(0,0) = vcparams(0);
		D(0,1) = D(1,0) = vcparams(1);
		D(1,1) = vcparams(2);
	/**** update D *****************************************************************************************************************************/

	/**** update V *****************************************************************************************************************************/
		for (int i=0; i<n; ++i)
		{
			V[i].setIdentity();
			V[i] *= vcparams(3);
			V[i].noalias() 	+= U.block(index_obs(i,0),0,index_obs(i,1),2)*D*U.block(index_obs(i,0),0,index_obs(i,1),2).transpose();
			V_inv[i] = V[i].selfadjointView<Eigen::Upper>().ldlt().solve(MatrixXd::Identity(index_obs(i,1),index_obs(i,1)));
		}
	/**** update V *****************************************************************************************************************************/

	/**** update P_theta ***********************************************************************************************************************/
		for (int i=0; i<n_minus_n2; ++i)
		{
			idx = i+n2;
			idxx = index_obs(idx,0)-nobs2;
			for (int k=0; k<m; ++k)
			{
				P_theta(i,k) = (resi.block(idxx,k,index_obs(idx,1),1).transpose()*V_inv[idx]*resi.block(idxx,k,index_obs(idx,1),1))(0,0);
			}
		}
		P_theta /= -2.;
		P_theta = P_theta.array().exp();
	/**** update P_theta ***********************************************************************************************************************/

	/**** update q, q_row_sum ******************************************************************************************************************/
		// for (int i=0; i<n_minus_n2; ++i)
		// {
		// 	for (int k=0; k<m; ++k)
		// 	{
		// 		q(i,k) = P_theta(i,k)*pB(Bspline_uni_ind(i+n2),k);
		// 	}
		// }
		for (int i = 0; i < n_minus_n2; ++i)
		{
			Bspline_Mat.row(i) = pB.row(Bspline_uni_ind(i+n2));
		}
		q = P_theta.array() * Bspline_Mat.array();
		q_row_sum = q.rowwise().sum();
		for (int i=0; i<n_minus_n2; ++i)
		{
			q.row(i) /= q_row_sum(i);
		}
	/**** update q, q_row_sum ******************************************************************************************************************/

	/**** E-step *******************************************************************************************************************************/


	/**** M-step *******************************************************************************************************************************/

	/**** update theta *************************************************************************************************************************/
		XinvVX.setZero();
		XinvVY.setZero();

		for (int i=0; i<n2; ++i)
		{
			V_inv_1[i] 	= V_inv[i].rowwise().sum();
			V_inv_T[i] 	= V_inv[i]*T.segment(index_obs(i,0),index_obs(i,1));
			V_inv_11(i) = V_inv_1[i].sum();
			V_inv_1T(i) = V_inv_T[i].sum();
			V_inv_TT(i) = V_inv_T[i].dot(T.segment(index_obs(i,0),index_obs(i,1)));
			V_inv_1Y 	= V_inv_1[i].dot(Y.segment(index_obs(i,0),index_obs(i,1)));
			V_inv_TY 	= V_inv_T[i].dot(Y.segment(index_obs(i,0),index_obs(i,1)));

			XinvVX(0,0) += V_inv_11(i);
			XinvVX.block(0,1,1,X_nc).noalias() 			+= V_inv_11(i)*X.row(i);
			XinvVX.block(0,nZ,1,Z_nc).noalias() 		+= V_inv_11(i)*Z.row(i);
			XinvVX(0,nT) += V_inv_1T(i);
			XinvVX.block(0,nXT,1,X_nc).noalias() 		+= V_inv_1T(i)*X.row(i);
			XinvVX.block(1,1,X_nc,X_nc).noalias() 		+= V_inv_11(i)*X.row(i).transpose()*X.row(i);
			XinvVX.block(1,nZ,X_nc,Z_nc).noalias() 		+= V_inv_11(i)*X.row(i).transpose()*Z.row(i);
			XinvVX.block(1,nT,X_nc,1).noalias() 		+= V_inv_1T(i)*X.row(i).transpose();
			XinvVX.block(1,nXT,X_nc,X_nc).noalias() 	+= V_inv_1T(i)*X.row(i).transpose()*X.row(i);
			XinvVX.block(nZ,nZ,Z_nc,Z_nc).noalias() 	+= V_inv_11(i)*Z.row(i).transpose()*Z.row(i);
			XinvVX.block(nZ,nT,Z_nc,1).noalias() 		+= V_inv_1T(i)*Z.row(i).transpose();
			XinvVX.block(nZ,nXT,Z_nc,X_nc).noalias() 	+= V_inv_1T(i)*Z.row(i).transpose()*X.row(i);
			XinvVX(nT,nT) += V_inv_TT(i);
			XinvVX.block(nT,nXT,1,X_nc).noalias() 		+= V_inv_TT(i)*X.row(i);
			XinvVX.block(nXT,nXT,X_nc,X_nc).noalias() 	+= V_inv_TT(i)*X.row(i).transpose()*X.row(i);

			XinvVY(0) += V_inv_1Y;
			XinvVY.segment(1,X_nc).noalias() 	+= V_inv_1Y*X.row(i).transpose();
			XinvVY.segment(nZ,Z_nc).noalias() 	+= V_inv_1Y*Z.row(i).transpose();
			XinvVY(nT) += V_inv_TY;
			XinvVY.segment(nXT,X_nc).noalias() 	+= V_inv_TY*X.row(i).transpose();

			if (ZT_flag)
			{
				XinvVX.block(0,nZT,1,ZT_nc).noalias() 		+= V_inv_1T(i)*ZT.row(i);
				XinvVX.block(1,nZT,X_nc,ZT_nc).noalias() 	+= V_inv_1T(i)*X.row(i).transpose()*ZT.row(i);
				XinvVX.block(nZ,nZT,Z_nc,ZT_nc).noalias() 	+= V_inv_1T(i)*Z.row(i).transpose()*ZT.row(i);
				XinvVX.block(nT,nZT,1,ZT_nc).noalias() 		+= V_inv_TT(i)*ZT.row(i);
				XinvVX.block(nXT,nZT,X_nc,ZT_nc).noalias() 	+= V_inv_TT(i)*X.row(i).transpose()*ZT.row(i);
				XinvVX.block(nZT,nZT,ZT_nc,ZT_nc).noalias() += V_inv_TT(i)*ZT.row(i).transpose()*ZT.row(i);

				XinvVY.segment(nZT,ZT_nc).noalias() += V_inv_TY*ZT.row(i).transpose();
			}
		}

		for (int i=0; i<n_minus_n2; ++i)
		{
			idx = i+n2;

			V_inv_1[idx] = V_inv[idx].rowwise().sum();
			V_inv_T[idx] = V_inv[idx]*T.segment(index_obs(idx,0),index_obs(idx,1));
			V_inv_11(idx) = V_inv_1[idx].sum();
			V_inv_1T(idx) = V_inv_T[idx].sum();
			V_inv_TT(idx) = V_inv_T[idx].dot(T.segment(index_obs(idx,0),index_obs(idx,1)));
			V_inv_1Y = V_inv_1[idx].dot(Y.segment(index_obs(idx,0),index_obs(idx,1)));
			V_inv_TY = V_inv_T[idx].dot(Y.segment(index_obs(idx,0),index_obs(idx,1)));

			XinvVX(0,0) 	+= V_inv_11(idx);
			XinvVX.block(0,nZ,1,Z_nc).noalias() 	+= V_inv_11(idx)*Z.row(idx);
			XinvVX(0,nT) 	+= V_inv_1T(idx);
			XinvVX.block(nZ,nZ,Z_nc,Z_nc).noalias() 	+= V_inv_11(idx)*Z.row(idx).transpose()*Z.row(idx);
			XinvVX.block(nZ,nT,Z_nc,1).noalias() 	+= V_inv_1T(idx)*Z.row(idx).transpose();
			XinvVX(nT,nT) 	+= V_inv_TT(idx);

			XinvVY(0) 	+= V_inv_1Y;
			XinvVY.segment(nZ,Z_nc).noalias() 	+= V_inv_1Y*Z.row(idx).transpose();
			XinvVY(nT) 	+= V_inv_TY;

			if (ZT_flag)
			{
				XinvVX.block(0,nZT,1,ZT_nc).noalias() 	+= V_inv_1T(idx)*ZT.row(idx);
				XinvVX.block(nZ,nZT,Z_nc,ZT_nc).noalias() 	+= V_inv_1T(idx)*Z.row(idx).transpose()*ZT.row(idx);
				XinvVX.block(nT,nZT,1,ZT_nc).noalias() 	+= V_inv_TT(idx)*ZT.row(idx);
				XinvVX.block(nZT,nZT,ZT_nc,ZT_nc).noalias() 	+= V_inv_TT(idx)*ZT.row(idx).transpose()*ZT.row(idx);

				XinvVY.segment(nZT,ZT_nc).noalias() 	+= V_inv_TY*ZT.row(idx).transpose();
			}

			for (int k=0; k<m; ++k)
			{
				XinvVX.block(0,1,1,X_nc).noalias() 	+= q(i,k)*V_inv_11(idx)*X_uni.row(k);
				XinvVX.block(0,nXT,1,X_nc).noalias() 	+= q(i,k)*V_inv_1T(idx)*X_uni.row(k);
				XinvVX.block(1,1,X_nc,X_nc).noalias() 	+= q(i,k)*V_inv_11(idx)*X_uni.row(k).transpose()*X_uni.row(k);
				XinvVX.block(1,nZ,X_nc,Z_nc).noalias() 	+= q(i,k)*V_inv_11(idx)*X_uni.row(k).transpose()*Z.row(idx);
				XinvVX.block(1,nT,X_nc,1).noalias() 	+= q(i,k)*V_inv_1T(idx)*X_uni.row(k).transpose();
				XinvVX.block(1,nXT,X_nc,X_nc).noalias() 	+= q(i,k)*V_inv_1T(idx)*X_uni.row(k).transpose()*X_uni.row(k);
				XinvVX.block(nZ,nXT,Z_nc,X_nc).noalias() 	+= q(i,k)*V_inv_1T(idx)*Z.row(idx).transpose()*X_uni.row(k);
				XinvVX.block(nT,nXT,1,X_nc).noalias() 	+= q(i,k)*V_inv_TT(idx)*X_uni.row(k);
				XinvVX.block(nXT,nXT,X_nc,X_nc).noalias() 	+= q(i,k)*V_inv_TT(idx)*X_uni.row(k).transpose()*X_uni.row(k);

				XinvVY.segment(1,X_nc).noalias() 	+= q(i,k)*V_inv_1Y*X_uni.row(k).transpose();
				XinvVY.segment(nXT,X_nc).noalias() 	+= q(i,k)*V_inv_TY*X_uni.row(k).transpose();

				if (ZT_flag)
				{
					XinvVX.block(1,nZT,X_nc,ZT_nc).noalias() 	+= q(i,k)*V_inv_1T(idx)*X_uni.row(k).transpose()*ZT.row(idx);
					XinvVX.block(nXT,nZT,X_nc,ZT_nc).noalias() 	+= q(i,k)*V_inv_TT(idx)*X_uni.row(k).transpose()*ZT.row(idx);
				}
			}
		}

		theta = XinvVX.selfadjointView<Eigen::Upper>().ldlt().solve(XinvVY);
	/**** update theta *************************************************************************************************************************/

	/**** update residuals *********************************************************************************************************************/
		X_uni_theta = X_uni*theta.segment(1,X_nc);
		Z_theta = Z*theta.segment(nZ,Z_nc);
		TX_uni_theta = X_uni*theta.segment(nXT,X_nc);
		if (ZT_flag)
		{
			TZ_theta = ZT*theta.tail(ZT_nc);
		}
		for (int i=0; i<n2; ++i)
		{
			resi2.segment(index_obs(i,0),index_obs(i,1)) = Y.segment(index_obs(i,0),index_obs(i,1));
			resi2.segment(index_obs(i,0),index_obs(i,1)).array() -= theta(0)+X_uni_theta(X_uni_ind(i))+Z_theta(i);
			resi2.segment(index_obs(i,0),index_obs(i,1)).noalias() -= T.segment(index_obs(i,0),index_obs(i,1))*theta(nT);
			resi2.segment(index_obs(i,0),index_obs(i,1)).noalias() -= T.segment(index_obs(i,0),index_obs(i,1))*TX_uni_theta(X_uni_ind(i));
			if (ZT_flag)
			{
				resi2.segment(index_obs(i,0),index_obs(i,1)).noalias() -= T.segment(index_obs(i,0),index_obs(i,1))*TZ_theta(i);
			}
		}
		for (int i=0; i<n_minus_n2; ++i)
		{
			idx = i+n2;
			idxx = index_obs(idx,0)-nobs2;
			resi.block(idxx,0,index_obs(idx,1),1) = Y.segment(index_obs(idx,0),index_obs(idx,1));
			resi.block(idxx,0,index_obs(idx,1),1).array() -= theta(0)+Z_theta(idx);
			resi.block(idxx,0,index_obs(idx,1),1).noalias() -= T.segment(index_obs(idx,0),index_obs(idx,1))*theta(nT);
			if (ZT_flag)
			{
				resi.block(idxx,0,index_obs(idx,1),1).noalias() -= T.segment(index_obs(idx,0),index_obs(idx,1))*TZ_theta(idx);
			}
			for (int k=1; k<m; ++k)
			{
				resi.block(idxx,k,index_obs(idx,1),1) = resi.block(idxx,0,index_obs(idx,1),1);
			}
			for (int k=0; k<m; ++k)
			{
				resi.block(idxx,k,index_obs(idx,1),1).array() -= X_uni_theta(k);
				resi.block(idxx,k,index_obs(idx,1),1).noalias() -= T.segment(index_obs(idx,0),index_obs(idx,1))*TX_uni_theta(k);
			}
		}
	/**** update residuals *********************************************************************************************************************/

	/**** update variance components parameters ************************************************************************************************/
		L_var.setZero();
		I_var.setZero();

		for (int i=0; i<n2; ++i)
		{
			V_inv_1R = V_inv_1[i].dot(resi2.segment(index_obs(i,0),index_obs(i,1)));
			V_inv_TR = V_inv_T[i].dot(resi2.segment(index_obs(i,0),index_obs(i,1)));
			V_inv_R2 = (V_inv[i]*resi2.segment(index_obs(i,0),index_obs(i,1))).squaredNorm();

			L_var(0) 	+= V_inv_1R*V_inv_1R-V_inv_11(i);
			L_var(1) 	+= 2.*(V_inv_TR*V_inv_1R-V_inv_1T(i));
			L_var(2) 	+= V_inv_TR*V_inv_TR-V_inv_TT(i);
			L_var(3) 	+= V_inv_R2-V_inv[i].trace();

			I_var(0,0) 	+= V_inv_11(i)*V_inv_11(i);
			I_var(0,1) 	+= 2.*V_inv_1T(i)*V_inv_11(i);
			I_var(0,2) 	+= V_inv_1T(i)*V_inv_1T(i);
			I_var(0,3) 	+= V_inv_1[i].squaredNorm();
			I_var(1,1) 	+= 2*(V_inv_1T(i)*V_inv_1T(i)+V_inv_11(i)*V_inv_TT(i));
			I_var(1,2) 	+= 2*V_inv_1T(i)*V_inv_TT(i);
			I_var(1,3) 	+= 2*V_inv_1[i].dot(V_inv_T[i]);
			I_var(2,2) 	+= V_inv_TT(i)*V_inv_TT(i);
			I_var(2,3) 	+= V_inv_T[i].squaredNorm();
			I_var(3,3) 	+= (V_inv[i]*V_inv[i]).trace();
		}

		for (int i=0; i<n_minus_n2; ++i)
		{
			idx = i+n2;
			idxx = index_obs(idx,0)-nobs2;

			for (int k=0; k<m; ++k)
			{
				// sum(V_inv[idx] * resi.block)
				V_inv_1R = (V_inv_1[idx].transpose()*resi.block(idxx,k,index_obs(idx,1),1))(0,0);
				V_inv_TR = (V_inv_T[idx].transpose()*resi.block(idxx,k,index_obs(idx,1),1))(0,0);
				V_inv_R2 = (V_inv[idx]*resi.block(idxx,k,index_obs(idx,1),1)).squaredNorm();

				L_var(0) 	+= q(i,k)*V_inv_1R*V_inv_1R;
				L_var(1) 	+= 2.*q(i,k)*V_inv_TR*V_inv_1R;
				L_var(2) 	+= q(i,k)*V_inv_TR*V_inv_TR;
				L_var(3) 	+= q(i,k)*V_inv_R2;
			}

			L_var(0) -= V_inv_11(idx);
			L_var(1) -= 2.*V_inv_1T(idx);
			L_var(2) -= V_inv_TT(idx);
			L_var(3) -= V_inv[idx].trace();

			I_var(0,0) 	+= V_inv_11(idx)*V_inv_11(idx);
			I_var(0,1) 	+= 2.*V_inv_1T(idx)*V_inv_11(idx);
			I_var(0,2) 	+= V_inv_1T(idx)*V_inv_1T(idx);
			I_var(0,3) 	+= V_inv_1[idx].squaredNorm();
			I_var(1,1) 	+= 2*(V_inv_1T(idx)*V_inv_1T(idx)+V_inv_11(idx)*V_inv_TT(idx));
			I_var(1,2) 	+= 2*V_inv_1T(idx)*V_inv_TT(idx);
			I_var(1,3) 	+= 2*V_inv_1[idx].dot(V_inv_T[idx]);
			I_var(2,2) 	+= V_inv_TT(idx)*V_inv_TT(idx);
			I_var(2,3) 	+= V_inv_T[idx].squaredNorm();
			I_var(3,3) 	+= (V_inv[idx]*V_inv[idx]).trace();
		}

		vcparams = I_var.selfadjointView<Eigen::Upper>().ldlt().solve(L_var);
		vcparams += vcparams0;
	/**** update variance components parameters ************************************************************************************************/

	/**** update p *****************************************************************************************************************************/
		p.setZero();
		pthetaOverQ = P_theta.array().colwise() / q_row_sum.array();

		for (int i=0; i<n_minus_n2; ++i)
		{
			// for (int k=0; k<m; ++k)
			// {
			// 	for (int j=0; j<s; ++j)
			// 	{
			// 		p(k,j) 	+= Bspline_uni(Bspline_uni_ind(i+n2),j)*P_theta(i,k)/q_row_sum(i);
			// 	}
			// }
			p += pthetaOverQ.row(i).transpose() * Bspline_uni.row(Bspline_uni_ind(i + n2));

		}
		p = p.array()*p0.array();
		p += p_static;
		p_col_sum = p.colwise().sum();
		for (int j=0; j<s; ++j)
		{
			p.col(j) /= p_col_sum(j);
		}
	/**** update p *****************************************************************************************************************************/

	/**** M-step *******************************************************************************************************************************/


	/**** calculate the sum of absolute differences between estimates in the current and previous iterations ***********************************/
		tol = (theta-theta0).array().abs().sum();
		tol += (vcparams-vcparams0).array().abs().sum();
		tol += (p-p0).array().abs().sum();
	/**** calculate the sum of absolute differences between estimates in the current and previous iterations ***********************************/

	/**** update parameters ********************************************************************************************************************/
		theta0 = theta;
		vcparams0 = vcparams;
		p0 = p;
	/**** update parameters ********************************************************************************************************************/

	/**** check convergence ********************************************************************************************************************/
		if (tol < TOL)
		{
			break;
		}
	/**** check convergence ********************************************************************************************************************/

	// /* RT's test code */
	// time(&t2);
	// Rcout << iter << '\t' << difftime(t2,t1) << '\t' << tol << endl;
	// /* RT's test code end */
	}
	/**** EM algorithm *****************************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** variance estimation **********************************************************************************************************************/
	if (iter == MAX_ITER)
	{
		flag_nonconvergence = true;
		flag_nonconvergence_cov = true;
		theta.setConstant(-999.);
		cov_theta.setConstant(-999.);
		vcparams.setConstant(-999.);
	}
	else if (noSE)
	{
		flag_nonconvergence_cov = true;
		cov_theta.setConstant(-999.);
	}
	else
	{
		VectorXd profile_vec(ncov);
		MatrixXd logp(m,s);
		MatrixXd profile_mat(ncov,ncov);
		double loglik;

		profile_mat.setZero();
		profile_vec.setZero();

		loglik = LMM_GeneralSplineProfile(pB,p_col_sum,q_row_sum,vcparams,vcparams0,p,p0,P_theta,q,logp,Z_theta,X_uni_theta,TX_uni_theta,resi2,resi,
			D,V,V_inv,V_inv_1,V_inv_T,V_inv_11,V_inv_1T,V_inv_TT,L_var,I_var,TZ_theta,theta,vc_initial,Y,T,U,index_obs,Bspline_uni,Z,X_uni,X_uni_ind,
			Bspline_uni_ind,p_static,n,n2,m,s,n_minus_n2,nobs,nobs2,X_nc,Z_nc,nT,MAX_ITER,TOL,nZ,nXT,ZT_flag,ZT_nc,ZT);
		if (loglik == -999.)
		{
			flag_nonconvergence_cov = true;
		}
		// for (int i=0; i<ncov; ++i)
		// {
		// 	for (int j=i; j<ncov; ++j)
		// 	{
		// 		profile_mat(i,j) = loglik;
		// 	}
		// }
		profile_mat.triangularView<Upper>().setConstant(loglik);

		for (int i=0; i<ncov; ++i)
		{
			theta0 = theta;
			theta0(i) += hn;
			profile_vec(i) = LMM_GeneralSplineProfile(pB,p_col_sum,q_row_sum,vcparams,vcparams0,p,p0,P_theta,q,logp,Z_theta,X_uni_theta,TX_uni_theta,resi2,resi,
				D,V,V_inv,V_inv_1,V_inv_T,V_inv_11,V_inv_1T,V_inv_TT,L_var,I_var,TZ_theta,theta0,vc_initial,Y,T,U,index_obs,Bspline_uni,Z,X_uni,X_uni_ind,
				Bspline_uni_ind,p_static,n,n2,m,s,n_minus_n2,nobs,nobs2,X_nc,Z_nc,nT,MAX_ITER,TOL,nZ,nXT,ZT_flag,ZT_nc,ZT);
			if(profile_vec(i) == -999.)
			{
				flag_nonconvergence_cov = true;
			}
		}

		for (int i=0; i<ncov; ++i)
		{
			for (int j=i; j<ncov; ++j)
			{
				theta0 = theta;
				theta0(i) += hn;
				theta0(j) += hn;
				loglik = LMM_GeneralSplineProfile(pB,p_col_sum,q_row_sum,vcparams,vcparams0,p,p0,P_theta,q,logp,Z_theta,X_uni_theta,TX_uni_theta,resi2,resi,
					D,V,V_inv,V_inv_1,V_inv_T,V_inv_11,V_inv_1T,V_inv_TT,L_var,I_var,TZ_theta,theta0,vc_initial,Y,T,U,index_obs,Bspline_uni,Z,X_uni,X_uni_ind,
					Bspline_uni_ind,p_static,n,n2,m,s,n_minus_n2,nobs,nobs2,X_nc,Z_nc,nT,MAX_ITER,TOL,nZ,nXT,ZT_flag,ZT_nc,ZT);
				if (loglik == -999.)
				{
					flag_nonconvergence_cov = true;
				}
				profile_mat(i,j) += loglik;
				profile_mat(i,j) -= profile_vec(i)+profile_vec(j);
			}
		}

		if (flag_nonconvergence_cov == true)
		{
			cov_theta.setConstant(-999.);
		}
		else
		{
			// for (int i=0; i<ncov; ++i)
			// {
			// 	for (int j=i+1; j<ncov; ++j)
			// 	{
			// 		profile_mat(j,i) = profile_mat(i,j);
			// 	}
			// }
			profile_mat = profile_mat.selfadjointView<Upper>();
			profile_mat /= hn*hn;
			profile_mat = -profile_mat;
			cov_theta = profile_mat.ldlt().solve(MatrixXd::Identity(ncov,ncov));
		}
	}
	/**** variance estimation **********************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** return output to R ***********************************************************************************************************************/
	return List::create(Named("theta") = theta,
		Named("cov_theta") = cov_theta,
		Named("vcparams") = vcparams,
		Named("flag_nonconvergence") = flag_nonconvergence,
		Named("flag_nonconvergence_cov") = flag_nonconvergence_cov);
	/**** return output to R ***********************************************************************************************************************/
	/*#############################################################################################################################################*/
} // LMM_GeneralSpline
