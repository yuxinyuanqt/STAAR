// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_STAAR_RCPPEXPORTS_H_GEN_
#define RCPP_STAAR_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace STAAR {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("STAAR", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("STAAR", "_STAAR_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in STAAR");
            }
        }
    }

    inline double Bisection(arma::vec egvalues, double q, double xmin, double xmax) {
        typedef SEXP(*Ptr_Bisection)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_Bisection p_Bisection = NULL;
        if (p_Bisection == NULL) {
            validateSignature("double(*Bisection)(arma::vec,double,double,double)");
            p_Bisection = (Ptr_Bisection)R_GetCCallable("STAAR", "_STAAR_Bisection");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_Bisection(Shield<SEXP>(Rcpp::wrap(egvalues)), Shield<SEXP>(Rcpp::wrap(q)), Shield<SEXP>(Rcpp::wrap(xmin)), Shield<SEXP>(Rcpp::wrap(xmax)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double CCT_pval(arma::vec x, arma::vec weights) {
        typedef SEXP(*Ptr_CCT_pval)(SEXP,SEXP);
        static Ptr_CCT_pval p_CCT_pval = NULL;
        if (p_CCT_pval == NULL) {
            validateSignature("double(*CCT_pval)(arma::vec,arma::vec)");
            p_CCT_pval = (Ptr_CCT_pval)R_GetCCallable("STAAR", "_STAAR_CCT_pval");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_CCT_pval(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(weights)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double K(double x, arma::vec egvalues) {
        typedef SEXP(*Ptr_K)(SEXP,SEXP);
        static Ptr_K p_K = NULL;
        if (p_K == NULL) {
            validateSignature("double(*K)(double,arma::vec)");
            p_K = (Ptr_K)R_GetCCallable("STAAR", "_STAAR_K");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_K(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(egvalues)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double K1(double x, arma::vec egvalues, double q) {
        typedef SEXP(*Ptr_K1)(SEXP,SEXP,SEXP);
        static Ptr_K1 p_K1 = NULL;
        if (p_K1 == NULL) {
            validateSignature("double(*K1)(double,arma::vec,double)");
            p_K1 = (Ptr_K1)R_GetCCallable("STAAR", "_STAAR_K1");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_K1(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(egvalues)), Shield<SEXP>(Rcpp::wrap(q)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double K2(double x, arma::vec egvalues) {
        typedef SEXP(*Ptr_K2)(SEXP,SEXP);
        static Ptr_K2 p_K2 = NULL;
        if (p_K2 == NULL) {
            validateSignature("double(*K2)(double,arma::vec)");
            p_K2 = (Ptr_K2)R_GetCCallable("STAAR", "_STAAR_K2");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_K2(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(egvalues)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double Saddle(double q, arma::vec egvalues) {
        typedef SEXP(*Ptr_Saddle)(SEXP,SEXP);
        static Ptr_Saddle p_Saddle = NULL;
        if (p_Saddle == NULL) {
            validateSignature("double(*Saddle)(double,arma::vec)");
            p_Saddle = (Ptr_Saddle)R_GetCCallable("STAAR", "_STAAR_Saddle");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_Saddle(Shield<SEXP>(Rcpp::wrap(q)), Shield<SEXP>(Rcpp::wrap(egvalues)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline List matrix_flip(arma::mat G) {
        typedef SEXP(*Ptr_matrix_flip)(SEXP);
        static Ptr_matrix_flip p_matrix_flip = NULL;
        if (p_matrix_flip == NULL) {
            validateSignature("List(*matrix_flip)(arma::mat)");
            p_matrix_flip = (Ptr_matrix_flip)R_GetCCallable("STAAR", "_STAAR_matrix_flip");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_matrix_flip(Shield<SEXP>(Rcpp::wrap(G)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline arma::mat matrix_impute(arma::mat G) {
        typedef SEXP(*Ptr_matrix_impute)(SEXP);
        static Ptr_matrix_impute p_matrix_impute = NULL;
        if (p_matrix_impute == NULL) {
            validateSignature("arma::mat(*matrix_impute)(arma::mat)");
            p_matrix_impute = (Ptr_matrix_impute)R_GetCCallable("STAAR", "_STAAR_matrix_impute");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_matrix_impute(Shield<SEXP>(Rcpp::wrap(G)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

}

#endif // RCPP_STAAR_RCPPEXPORTS_H_GEN_
