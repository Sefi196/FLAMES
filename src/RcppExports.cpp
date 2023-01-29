// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// find_barcode
Rcpp::List find_barcode(Rcpp::String fastq_dir, Rcpp::String stats_file, Rcpp::String out_fastq, Rcpp::String ref_csv, int MAX_DIST, int UMI_LEN);
RcppExport SEXP _FLAMES_find_barcode(SEXP fastq_dirSEXP, SEXP stats_fileSEXP, SEXP out_fastqSEXP, SEXP ref_csvSEXP, SEXP MAX_DISTSEXP, SEXP UMI_LENSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type fastq_dir(fastq_dirSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type stats_file(stats_fileSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type out_fastq(out_fastqSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type ref_csv(ref_csvSEXP);
    Rcpp::traits::input_parameter< int >::type MAX_DIST(MAX_DISTSEXP);
    Rcpp::traits::input_parameter< int >::type UMI_LEN(UMI_LENSEXP);
    rcpp_result_gen = Rcpp::wrap(find_barcode(fastq_dir, stats_file, out_fastq, ref_csv, MAX_DIST, UMI_LEN));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FLAMES_find_barcode", (DL_FUNC) &_FLAMES_find_barcode, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_FLAMES(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
