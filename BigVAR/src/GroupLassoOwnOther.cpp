#include <Rcpp.h>
#include <RcppEigen.h>

Eigen::MatrixXd getRows(const Eigen::MatrixXd &mat,
                        const Rcpp::IntegerVector &rows) {
    Eigen::MatrixXd ret(rows.size(), mat.cols());
    for (int i = 0; i < rows.size(); i++) {
        ret.row(i) = mat.row(rows(i));
    }
    return(ret);
}

Eigen::MatrixXd getCols(const Eigen::MatrixXd &mat,
                        const Rcpp::IntegerVector &cols) {
    Eigen::MatrixXd ret(mat.rows(), cols.size());
    for (int i = 0; i < cols.size(); i++) {
        ret.col(i) = mat.col(cols(i));
    }
    return(ret);
}

Eigen::VectorXd matvecProd(const Eigen::MatrixXd &mat,
                           const Eigen::VectorXd &vec,
                           const Rcpp::IntegerVector &idx) {
    Eigen::VectorXd ret = Eigen::VectorXd::Zero(mat.rows());
    for (int i = 0; i < idx.size(); i++) {
        ret.array() += mat.col(idx(i)).array() * vec(idx(i));
    }
    return(ret);
}

Eigen::VectorXd getElems(const Eigen::VectorXd &vec,
                         const Rcpp::IntegerVector &elems) {
    Eigen::VectorXd ret(elems.size());
    for (int i = 0; i < elems.size(); i++) {
        ret(i) = vec(elems(i));
    }
    return(ret);
}

void setElements(Eigen::VectorXd &vec,
                            const Rcpp::IntegerVector &idx,
                            const Eigen::VectorXd &elems) {
    for (int i = 0; i < idx.size(); i++) {
        vec(idx(i)) = elems(i);
    }
}

// Newton-Raphson Functions
double trust32 (int k, const Eigen::MatrixXd& PP, double delta, double lambda, const Eigen::VectorXd& EigVA, const Eigen::MatrixXd& EigVector) {
	double g = 0;
	for (int i = 0; i < k; i++) {
        g += std::pow((EigVector.col(i).transpose() * PP)(0), 2) / std::pow((EigVA(i) * delta + lambda), 2);
    }
	return(g);
}

double fprime2 (int k, const Eigen::MatrixXd& PP, double delta, double lambda, const Eigen::VectorXd& EigVA, const Eigen::MatrixXd& EigVE) {
	double gg2 = 0;
	for (int i = 0; i < k; ++i) {
        gg2 += (std::pow((EigVE.col(i).transpose() * PP)(0), 2) * EigVA(i)) / std::pow((EigVA(i) * delta + lambda), 3);
    }
	double c1 = trust32(k, PP, delta, lambda, EigVA, EigVE);
	double res = -.5 * std::pow(c1, -1.5) * -2 * gg2;
	return(res);
}

//Newton raphson for trust region problem
double Newton2 (int k, const Eigen::MatrixXd& P, double lambda, const Eigen::VectorXd& EigVA, const Eigen::MatrixXd& EigVE){
	double delta = 0;
	double threshold = 1;
	double phi = 0;
	double deltanew = delta;
	int iter = 0;
	while (threshold > 0.0001 && iter < 20) {
        phi = 1 - 1 / std::pow(trust32(k, P, delta, lambda, EigVA, EigVE), 0.5);
        deltanew += phi / fprime2(k, P, deltanew, lambda, EigVA, EigVE);
        threshold = fabs(delta - deltanew);
        delta = deltanew;
        phi = 0;
        iter += 1;
    }
	return(deltanew);
}

Rcpp::List BlockUpdate2(const Eigen::MatrixXd& ZZ1, double lam, const Eigen::VectorXd& Y1, double eps, Rcpp::List groups, Rcpp::List fullgroups, Rcpp::List compgroups, int k, Rcpp::List M2f_, Rcpp::List eigvalF_, Rcpp::List eigvecF_, Eigen::VectorXd& B, int k1, int& converge){
	int n1 = groups.size();
	Rcpp::List active(n1);
	Eigen::VectorXd BPrev = B;
	int count = 0;

	if (groups.size() == count) {
        B.fill(0);
        active = groups;
	} else {
		for (int i=0; i < n1; ++i)	{
            Rcpp::IntegerVector s1 = groups[i];
            Rcpp::IntegerVector s2 = fullgroups[i];
            Rcpp::IntegerVector scomp = compgroups[i];

            if (max(s1) == 0) {
                setElements(B, s2, Eigen::VectorXd::Zero(s2.size()));
                active(i) = 0;
            }
            if (max(s1) != 0) {
                Eigen::VectorXd r = matvecProd(ZZ1, B, scomp) - Y1;

                Eigen::MatrixXd M1 = getCols(ZZ1, s1);
                Eigen::MatrixXd M2 = M2f_(i);
                Eigen::VectorXd eigval = eigvalF_(i);
                Eigen::MatrixXd eigvec = eigvecF_(i);
                Eigen::MatrixXd p = M1.transpose() * r;

                double rho = sqrt(static_cast<double>(s1.size()));
                double adjlam = rho*lam;

                if (p.lpNorm<2>() <= adjlam) {
                    active(i) = 0;
                } else {
                    int k1a = M2.cols();
                    double deltfin = Newton2(k1a, p, adjlam, eigval, eigvec);

                    Eigen::MatrixXd D1 = Eigen::MatrixXd::Identity(s1.size(), s1.size());

                    //correct for rare occurrence where newton returns zero
                    if (deltfin == 0) {
                        deltfin += std::numeric_limits<double>::epsilon();
                    }

                    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PQR;
                    PQR.compute(M2.array() + adjlam / deltfin * D1.array());

                    Eigen::VectorXd astar = -PQR.solve(p);
                    setElements(B, s1, astar);
                    active(i) = s1;
                }
            }
        }
	}

	double thresh = (B - BPrev).lpNorm<Eigen::Infinity>();
	if (thresh < eps) {
        converge = 1;
    } else {
		converge = 0;
	}

	return(active);
}

Eigen::VectorXd ThreshUpdateOO(const Eigen::MatrixXd& ZZ,
                               double lam,
                               const Eigen::VectorXd& Y,
                               double eps,
                               Rcpp::List groups,
                               Rcpp::List fullgroups,
                               Rcpp::List compgroups,
                               Rcpp::List M2f_,
                               Rcpp::List eigvalF_,
                               Rcpp::List eigvecF_,
                               Eigen::VectorXd& B,
                               int n,
                               int k1)
{
    int converge = 0;
    int n1 = groups.size();
    Eigen::VectorXd BPrev = B;
    Rcpp::List active(n1);
    int count=0;
    for (int i=0; i < n1; ++i) {
        Rcpp::NumericVector g1 = groups[i];
        count += max(g1);
    }
    if (count == 0) {
        B.fill(0.0);
        active = groups;
    } else {
        double threshold = 10 * eps;
        while (threshold > eps) {
            active = BlockUpdate2(ZZ, lam, Y, eps, groups, fullgroups, compgroups, n, M2f_, eigvalF_, eigvecF_, B, k1, converge);
            threshold = (B - BPrev).lpNorm<Eigen::Infinity>();
            BPrev = B;
      }
    }
    return(B);
}

//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::export]]
Rcpp::List GamLoopGLOO(Rcpp::NumericVector beta_,
                       Rcpp::List Activeset,
                       Rcpp::NumericVector gamm,
                       const Eigen::VectorXd& Y,
                       const Eigen::MatrixXd& Z,
                       Rcpp::List jj,
                       Rcpp::List jjfull,
                       Rcpp::List jjcomp,
                       double eps,
                       Eigen::VectorXd& YMean2,
                       Eigen::VectorXd& ZMean2,
                       int k,
                       int pk,
                       Rcpp::List M2f_,
                       Rcpp::List eigvalF_,
                       Rcpp::List eigvecF_,
                       int k1)
{
	Rcpp::IntegerVector dims = beta_.attr("dim");
	int gran2 = gamm.size();
	Rcpp::List activefinal(gran2);
	Eigen::Map<Eigen::MatrixXd> beta2 = Eigen::Map<Eigen::MatrixXd>(beta_.begin(), dims[2], dims[1]*dims[0]);
	Eigen::MatrixXd betafin = Eigen::MatrixXd::Zero(dims[0]*dims[2], dims[1]+1);

	Rcpp::List iterations(gran2);

    Eigen::VectorXd B = Eigen::VectorXd(dims[0]*dims[1]);

	for (int i=0; i < gran2; ++i) {
        double gam = gamm[i];
        B = beta2.row(i);
        Rcpp::List Active = Activeset[i];
        int k2 = 0;
        int converge = 0;
        while (converge == 0) {
            B = ThreshUpdateOO(Z, gam, Y, eps, Active, jjfull, jjcomp, M2f_, eigvalF_, eigvecF_, B, k1, k1);
            Active = BlockUpdate2(Z, gam, Y, eps, jjfull, jjfull, jjcomp, k, M2f_, eigvalF_, eigvecF_, B, k1, converge);
            k2++;
        }
        Eigen::MatrixXd betaF = Eigen::Map<Eigen::MatrixXd>(B.data(), dims[0], dims[1]);
        Eigen::VectorXd nu = YMean2 - betaF * ZMean2;
        Eigen::MatrixXd betafin_slice = Eigen::MatrixXd(nu.rows(), nu.cols()+betaF.cols());
        betafin_slice << nu, betaF;
        betafin.middleRows(i*dims[0], dims[0]) = betafin_slice;
        activefinal[i] = Active;
        iterations[i] = k2;
    }
	Rcpp::List Results = Rcpp::List::create(Rcpp::Named("beta")=betafin, Rcpp::Named("active")=wrap(activefinal), Rcpp::Named("iterations")=iterations);
	return(Results);
}
