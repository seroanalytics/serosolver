// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// antibody_model
NumericVector antibody_model(const NumericMatrix theta, const IntegerVector& unique_theta_indices, const IntegerVector& unique_biomarker_groups, const IntegerMatrix& infection_history_mat, const IntegerVector& infection_history_mat_indices, const IntegerVector& indiv_theta_groups, const NumericVector& possible_exposure_times, const IntegerVector& possible_exposure_times_indices, const IntegerVector& exposure_groups, const IntegerVector& unique_exposure_groups, const NumericVector& sample_times, const IntegerVector& type_data_start, const IntegerVector& biomarker_groups, const IntegerVector& sample_data_start, const IntegerVector& antibody_data_start, const IntegerVector& nrows_per_sample, const IntegerVector& biomarker_id_indices, const IntegerVector& start_level_indices, const NumericVector& starting_antibody_levels, const NumericVector& births, const arma::cube& antigenic_map_long, const arma::cube& antigenic_map_short, const NumericVector& antigenic_distances, const bool timevarying_groups, const bool variant_specific_pars, bool boost_before_infection);
RcppExport SEXP _serosolver_antibody_model(SEXP thetaSEXP, SEXP unique_theta_indicesSEXP, SEXP unique_biomarker_groupsSEXP, SEXP infection_history_matSEXP, SEXP infection_history_mat_indicesSEXP, SEXP indiv_theta_groupsSEXP, SEXP possible_exposure_timesSEXP, SEXP possible_exposure_times_indicesSEXP, SEXP exposure_groupsSEXP, SEXP unique_exposure_groupsSEXP, SEXP sample_timesSEXP, SEXP type_data_startSEXP, SEXP biomarker_groupsSEXP, SEXP sample_data_startSEXP, SEXP antibody_data_startSEXP, SEXP nrows_per_sampleSEXP, SEXP biomarker_id_indicesSEXP, SEXP start_level_indicesSEXP, SEXP starting_antibody_levelsSEXP, SEXP birthsSEXP, SEXP antigenic_map_longSEXP, SEXP antigenic_map_shortSEXP, SEXP antigenic_distancesSEXP, SEXP timevarying_groupsSEXP, SEXP variant_specific_parsSEXP, SEXP boost_before_infectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const NumericMatrix >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type unique_theta_indices(unique_theta_indicesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type unique_biomarker_groups(unique_biomarker_groupsSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type infection_history_mat(infection_history_matSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type infection_history_mat_indices(infection_history_mat_indicesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type indiv_theta_groups(indiv_theta_groupsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type possible_exposure_times(possible_exposure_timesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type possible_exposure_times_indices(possible_exposure_times_indicesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type exposure_groups(exposure_groupsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type unique_exposure_groups(unique_exposure_groupsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sample_times(sample_timesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type type_data_start(type_data_startSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type biomarker_groups(biomarker_groupsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type sample_data_start(sample_data_startSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type antibody_data_start(antibody_data_startSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type nrows_per_sample(nrows_per_sampleSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type biomarker_id_indices(biomarker_id_indicesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type start_level_indices(start_level_indicesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type starting_antibody_levels(starting_antibody_levelsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type births(birthsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type antigenic_map_long(antigenic_map_longSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type antigenic_map_short(antigenic_map_shortSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type antigenic_distances(antigenic_distancesSEXP);
    Rcpp::traits::input_parameter< const bool >::type timevarying_groups(timevarying_groupsSEXP);
    Rcpp::traits::input_parameter< const bool >::type variant_specific_pars(variant_specific_parsSEXP);
    Rcpp::traits::input_parameter< bool >::type boost_before_infection(boost_before_infectionSEXP);
    rcpp_result_gen = Rcpp::wrap(antibody_model(theta, unique_theta_indices, unique_biomarker_groups, infection_history_mat, infection_history_mat_indices, indiv_theta_groups, possible_exposure_times, possible_exposure_times_indices, exposure_groups, unique_exposure_groups, sample_times, type_data_start, biomarker_groups, sample_data_start, antibody_data_start, nrows_per_sample, biomarker_id_indices, start_level_indices, starting_antibody_levels, births, antigenic_map_long, antigenic_map_short, antigenic_distances, timevarying_groups, variant_specific_pars, boost_before_infection));
    return rcpp_result_gen;
END_RCPP
}
// subset_nullable_vector
NumericVector subset_nullable_vector(const Nullable<NumericVector>& x, int index1, int index2);
RcppExport SEXP _serosolver_subset_nullable_vector(SEXP xSEXP, SEXP index1SEXP, SEXP index2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type index1(index1SEXP);
    Rcpp::traits::input_parameter< int >::type index2(index2SEXP);
    rcpp_result_gen = Rcpp::wrap(subset_nullable_vector(x, index1, index2));
    return rcpp_result_gen;
END_RCPP
}
// logistic_transform
double logistic_transform(double p);
RcppExport SEXP _serosolver_logistic_transform(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(logistic_transform(p));
    return rcpp_result_gen;
END_RCPP
}
// logit_transform
double logit_transform(double p);
RcppExport SEXP _serosolver_logit_transform(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(logit_transform(p));
    return rcpp_result_gen;
END_RCPP
}
// transform_parameters_cpp
NumericMatrix transform_parameters_cpp(NumericVector pars, List scale_table, IntegerVector theta_indices, IntegerVector scale_par_indices, NumericMatrix demographics, IntegerVector transforms);
RcppExport SEXP _serosolver_transform_parameters_cpp(SEXP parsSEXP, SEXP scale_tableSEXP, SEXP theta_indicesSEXP, SEXP scale_par_indicesSEXP, SEXP demographicsSEXP, SEXP transformsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< List >::type scale_table(scale_tableSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type theta_indices(theta_indicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type scale_par_indices(scale_par_indicesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type demographics(demographicsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type transforms(transformsSEXP);
    rcpp_result_gen = Rcpp::wrap(transform_parameters_cpp(pars, scale_table, theta_indices, scale_par_indices, demographics, transforms));
    return rcpp_result_gen;
END_RCPP
}
// get_starting_antibody_levels
NumericVector get_starting_antibody_levels(const int n_measurements, const double min_measurement, const Nullable<NumericVector>& starting_antibody_levels);
RcppExport SEXP _serosolver_get_starting_antibody_levels(SEXP n_measurementsSEXP, SEXP min_measurementSEXP, SEXP starting_antibody_levelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n_measurements(n_measurementsSEXP);
    Rcpp::traits::input_parameter< const double >::type min_measurement(min_measurementSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type starting_antibody_levels(starting_antibody_levelsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_starting_antibody_levels(n_measurements, min_measurement, starting_antibody_levels));
    return rcpp_result_gen;
END_RCPP
}
// sum_likelihoods
NumericVector sum_likelihoods(NumericVector liks, IntegerVector indices, int n_indivs);
RcppExport SEXP _serosolver_sum_likelihoods(SEXP liksSEXP, SEXP indicesSEXP, SEXP n_indivsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type liks(liksSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< int >::type n_indivs(n_indivsSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_likelihoods(liks, indices, n_indivs));
    return rcpp_result_gen;
END_RCPP
}
// create_cross_reactivity_vector
NumericVector create_cross_reactivity_vector(NumericVector x, double cr_gradient);
RcppExport SEXP _serosolver_create_cross_reactivity_vector(SEXP xSEXP, SEXP cr_gradientSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type cr_gradient(cr_gradientSEXP);
    rcpp_result_gen = Rcpp::wrap(create_cross_reactivity_vector(x, cr_gradient));
    return rcpp_result_gen;
END_RCPP
}
// sum_buckets
NumericVector sum_buckets(NumericVector a, NumericVector buckets);
RcppExport SEXP _serosolver_sum_buckets(SEXP aSEXP, SEXP bucketsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type buckets(bucketsSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_buckets(a, buckets));
    return rcpp_result_gen;
END_RCPP
}
// sum_infections_by_group
IntegerMatrix sum_infections_by_group(IntegerMatrix inf_hist, NumericVector group_ids_vec, int n_groups, bool timevarying_groups);
RcppExport SEXP _serosolver_sum_infections_by_group(SEXP inf_histSEXP, SEXP group_ids_vecSEXP, SEXP n_groupsSEXP, SEXP timevarying_groupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type inf_hist(inf_histSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type group_ids_vec(group_ids_vecSEXP);
    Rcpp::traits::input_parameter< int >::type n_groups(n_groupsSEXP);
    Rcpp::traits::input_parameter< bool >::type timevarying_groups(timevarying_groupsSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_infections_by_group(inf_hist, group_ids_vec, n_groups, timevarying_groups));
    return rcpp_result_gen;
END_RCPP
}
// add_measurement_shifts
void add_measurement_shifts(NumericVector& predicted_antibody_levels, const NumericVector& to_add, const int& start_index_in_data, const int& end_index_in_data);
RcppExport SEXP _serosolver_add_measurement_shifts(SEXP predicted_antibody_levelsSEXP, SEXP to_addSEXP, SEXP start_index_in_dataSEXP, SEXP end_index_in_dataSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type predicted_antibody_levels(predicted_antibody_levelsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type to_add(to_addSEXP);
    Rcpp::traits::input_parameter< const int& >::type start_index_in_data(start_index_in_dataSEXP);
    Rcpp::traits::input_parameter< const int& >::type end_index_in_data(end_index_in_dataSEXP);
    add_measurement_shifts(predicted_antibody_levels, to_add, start_index_in_data, end_index_in_data);
    return R_NilValue;
END_RCPP
}
// inf_mat_prior_cpp
double inf_mat_prior_cpp(const IntegerMatrix& infection_history, const IntegerVector& n_alive, double shape1, double shape2);
RcppExport SEXP _serosolver_inf_mat_prior_cpp(SEXP infection_historySEXP, SEXP n_aliveSEXP, SEXP shape1SEXP, SEXP shape2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type infection_history(infection_historySEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type n_alive(n_aliveSEXP);
    Rcpp::traits::input_parameter< double >::type shape1(shape1SEXP);
    Rcpp::traits::input_parameter< double >::type shape2(shape2SEXP);
    rcpp_result_gen = Rcpp::wrap(inf_mat_prior_cpp(infection_history, n_alive, shape1, shape2));
    return rcpp_result_gen;
END_RCPP
}
// inf_mat_prior_cpp_vector
double inf_mat_prior_cpp_vector(const IntegerMatrix& infection_history, const IntegerVector& n_alive, const NumericVector& shape1s, const NumericVector& shape2s);
RcppExport SEXP _serosolver_inf_mat_prior_cpp_vector(SEXP infection_historySEXP, SEXP n_aliveSEXP, SEXP shape1sSEXP, SEXP shape2sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type infection_history(infection_historySEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type n_alive(n_aliveSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type shape1s(shape1sSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type shape2s(shape2sSEXP);
    rcpp_result_gen = Rcpp::wrap(inf_mat_prior_cpp_vector(infection_history, n_alive, shape1s, shape2s));
    return rcpp_result_gen;
END_RCPP
}
// inf_mat_prior_group_cpp
double inf_mat_prior_group_cpp(const IntegerMatrix& n_infections, const IntegerMatrix& n_alive, const double shape1, const double shape2);
RcppExport SEXP _serosolver_inf_mat_prior_group_cpp(SEXP n_infectionsSEXP, SEXP n_aliveSEXP, SEXP shape1SEXP, SEXP shape2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type n_infections(n_infectionsSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type n_alive(n_aliveSEXP);
    Rcpp::traits::input_parameter< const double >::type shape1(shape1SEXP);
    Rcpp::traits::input_parameter< const double >::type shape2(shape2SEXP);
    rcpp_result_gen = Rcpp::wrap(inf_mat_prior_group_cpp(n_infections, n_alive, shape1, shape2));
    return rcpp_result_gen;
END_RCPP
}
// inf_mat_prior_group_cpp_vector
double inf_mat_prior_group_cpp_vector(const IntegerMatrix& n_infections, const IntegerMatrix& n_alive, const NumericVector& shape1s, const NumericVector& shape2s);
RcppExport SEXP _serosolver_inf_mat_prior_group_cpp_vector(SEXP n_infectionsSEXP, SEXP n_aliveSEXP, SEXP shape1sSEXP, SEXP shape2sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type n_infections(n_infectionsSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type n_alive(n_aliveSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type shape1s(shape1sSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type shape2s(shape2sSEXP);
    rcpp_result_gen = Rcpp::wrap(inf_mat_prior_group_cpp_vector(n_infections, n_alive, shape1s, shape2s));
    return rcpp_result_gen;
END_RCPP
}
// inf_mat_prior_total_group_cpp
double inf_mat_prior_total_group_cpp(const IntegerVector& n_infections_group, const IntegerVector& n_alive_group, double shape1, double shape2);
RcppExport SEXP _serosolver_inf_mat_prior_total_group_cpp(SEXP n_infections_groupSEXP, SEXP n_alive_groupSEXP, SEXP shape1SEXP, SEXP shape2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type n_infections_group(n_infections_groupSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type n_alive_group(n_alive_groupSEXP);
    Rcpp::traits::input_parameter< double >::type shape1(shape1SEXP);
    Rcpp::traits::input_parameter< double >::type shape2(shape2SEXP);
    rcpp_result_gen = Rcpp::wrap(inf_mat_prior_total_group_cpp(n_infections_group, n_alive_group, shape1, shape2));
    return rcpp_result_gen;
END_RCPP
}
// likelihood_func_fast
NumericVector likelihood_func_fast(const NumericVector& theta, const NumericVector& obs, const NumericVector& predicted_antibody_levels);
RcppExport SEXP _serosolver_likelihood_func_fast(SEXP thetaSEXP, SEXP obsSEXP, SEXP predicted_antibody_levelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type predicted_antibody_levels(predicted_antibody_levelsSEXP);
    rcpp_result_gen = Rcpp::wrap(likelihood_func_fast(theta, obs, predicted_antibody_levels));
    return rcpp_result_gen;
END_RCPP
}
// likelihood_func_fast_continuous
NumericVector likelihood_func_fast_continuous(const NumericVector& theta, const NumericVector& obs, const NumericVector& predicted_antibody_levels);
RcppExport SEXP _serosolver_likelihood_func_fast_continuous(SEXP thetaSEXP, SEXP obsSEXP, SEXP predicted_antibody_levelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type predicted_antibody_levels(predicted_antibody_levelsSEXP);
    rcpp_result_gen = Rcpp::wrap(likelihood_func_fast_continuous(theta, obs, predicted_antibody_levels));
    return rcpp_result_gen;
END_RCPP
}
// likelihood_func_fast_continuous_fp
NumericVector likelihood_func_fast_continuous_fp(const NumericVector& theta, const NumericVector& obs, const NumericVector& predicted_antibody_levels);
RcppExport SEXP _serosolver_likelihood_func_fast_continuous_fp(SEXP thetaSEXP, SEXP obsSEXP, SEXP predicted_antibody_levelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type predicted_antibody_levels(predicted_antibody_levelsSEXP);
    rcpp_result_gen = Rcpp::wrap(likelihood_func_fast_continuous_fp(theta, obs, predicted_antibody_levels));
    return rcpp_result_gen;
END_RCPP
}
// inf_hist_prop_prior_v3
arma::mat inf_hist_prop_prior_v3(arma::mat infection_history_mat, const IntegerVector& sampled_indivs, const IntegerVector& age_mask, const IntegerVector& sample_mask, const IntegerVector& proposal_inf_hist_distances, const IntegerVector& n_infs, double shape1, double shape2, const NumericVector& rand_ns, const double& proposal_inf_hist_indiv_swap_ratio);
RcppExport SEXP _serosolver_inf_hist_prop_prior_v3(SEXP infection_history_matSEXP, SEXP sampled_indivsSEXP, SEXP age_maskSEXP, SEXP sample_maskSEXP, SEXP proposal_inf_hist_distancesSEXP, SEXP n_infsSEXP, SEXP shape1SEXP, SEXP shape2SEXP, SEXP rand_nsSEXP, SEXP proposal_inf_hist_indiv_swap_ratioSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type infection_history_mat(infection_history_matSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type sampled_indivs(sampled_indivsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type age_mask(age_maskSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type sample_mask(sample_maskSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type proposal_inf_hist_distances(proposal_inf_hist_distancesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type n_infs(n_infsSEXP);
    Rcpp::traits::input_parameter< double >::type shape1(shape1SEXP);
    Rcpp::traits::input_parameter< double >::type shape2(shape2SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rand_ns(rand_nsSEXP);
    Rcpp::traits::input_parameter< const double& >::type proposal_inf_hist_indiv_swap_ratio(proposal_inf_hist_indiv_swap_ratioSEXP);
    rcpp_result_gen = Rcpp::wrap(inf_hist_prop_prior_v3(infection_history_mat, sampled_indivs, age_mask, sample_mask, proposal_inf_hist_distances, n_infs, shape1, shape2, rand_ns, proposal_inf_hist_indiv_swap_ratio));
    return rcpp_result_gen;
END_RCPP
}
// inf_hist_prop_prior_v2_and_v4
List inf_hist_prop_prior_v2_and_v4(const NumericMatrix& theta, const IntegerVector& unique_theta_indices, const IntegerVector& unique_biomarker_groups, const IntegerVector& indiv_group_indices, const IntegerMatrix& infection_history_mat, const IntegerVector& infection_history_mat_indices, const NumericVector& likelihoods_pre_proposal, const IntegerVector& sampled_indivs, const IntegerVector& n_times_samp_vec, const IntegerVector& age_mask, const IntegerVector& sample_mask, const IntegerMatrix& n_alive, IntegerMatrix& n_infections, IntegerVector& n_infected_group, const arma::cube& prior_lookup, const double& proposal_inf_hist_indiv_swap_ratio, const int& swap_distance, const bool& propose_from_prior, const double& shape1, const double& shape2, const NumericVector& possible_exposure_times, const IntegerVector& possible_exposure_times_indices, const IntegerVector& exposure_groups, const IntegerVector& unique_exposure_groups, const NumericVector& sample_times, const IntegerVector& type_data_start, const IntegerVector& biomarker_groups, const IntegerVector& sample_data_start, const IntegerVector& antibody_data_start, const IntegerVector& nrows_per_sample, const IntegerVector& cum_nrows_per_individual_in_data, const IntegerVector& cum_nrows_per_individual_in_repeat_data, const IntegerVector& popn_group_id_vec, const IntegerVector& biomarker_id_indices, const IntegerVector& start_level_indices, const NumericVector& starting_antibody_levels, const NumericVector& births, const arma::cube& antigenic_map_long, const arma::cube& antigenic_map_short, const NumericVector& antigenic_distances, const NumericVector& antibody_data, const NumericVector& antibody_data_repeats, const int& n_measurements_total, const IntegerVector& repeat_indices, const bool& repeat_data_exist, const NumericVector& measurement_shifts, IntegerVector proposal_iter, IntegerVector accepted_iter, IntegerVector proposal_swap, IntegerVector accepted_swap, IntegerMatrix overall_swap_proposals, IntegerMatrix overall_add_proposals, const NumericVector time_sample_probs, const IntegerVector& total_alive, const IntegerVector& data_types, const NumericVector& obs_weights, const IntegerVector& indiv_possible_exposure_times_indices, const IntegerVector& indiv_poss_exp_times_start, const IntegerVector& indiv_poss_exp_times_end, const bool timevarying_groups, const bool variant_specific_pars, const double temp, bool solve_likelihood);
RcppExport SEXP _serosolver_inf_hist_prop_prior_v2_and_v4(SEXP thetaSEXP, SEXP unique_theta_indicesSEXP, SEXP unique_biomarker_groupsSEXP, SEXP indiv_group_indicesSEXP, SEXP infection_history_matSEXP, SEXP infection_history_mat_indicesSEXP, SEXP likelihoods_pre_proposalSEXP, SEXP sampled_indivsSEXP, SEXP n_times_samp_vecSEXP, SEXP age_maskSEXP, SEXP sample_maskSEXP, SEXP n_aliveSEXP, SEXP n_infectionsSEXP, SEXP n_infected_groupSEXP, SEXP prior_lookupSEXP, SEXP proposal_inf_hist_indiv_swap_ratioSEXP, SEXP swap_distanceSEXP, SEXP propose_from_priorSEXP, SEXP shape1SEXP, SEXP shape2SEXP, SEXP possible_exposure_timesSEXP, SEXP possible_exposure_times_indicesSEXP, SEXP exposure_groupsSEXP, SEXP unique_exposure_groupsSEXP, SEXP sample_timesSEXP, SEXP type_data_startSEXP, SEXP biomarker_groupsSEXP, SEXP sample_data_startSEXP, SEXP antibody_data_startSEXP, SEXP nrows_per_sampleSEXP, SEXP cum_nrows_per_individual_in_dataSEXP, SEXP cum_nrows_per_individual_in_repeat_dataSEXP, SEXP popn_group_id_vecSEXP, SEXP biomarker_id_indicesSEXP, SEXP start_level_indicesSEXP, SEXP starting_antibody_levelsSEXP, SEXP birthsSEXP, SEXP antigenic_map_longSEXP, SEXP antigenic_map_shortSEXP, SEXP antigenic_distancesSEXP, SEXP antibody_dataSEXP, SEXP antibody_data_repeatsSEXP, SEXP n_measurements_totalSEXP, SEXP repeat_indicesSEXP, SEXP repeat_data_existSEXP, SEXP measurement_shiftsSEXP, SEXP proposal_iterSEXP, SEXP accepted_iterSEXP, SEXP proposal_swapSEXP, SEXP accepted_swapSEXP, SEXP overall_swap_proposalsSEXP, SEXP overall_add_proposalsSEXP, SEXP time_sample_probsSEXP, SEXP total_aliveSEXP, SEXP data_typesSEXP, SEXP obs_weightsSEXP, SEXP indiv_possible_exposure_times_indicesSEXP, SEXP indiv_poss_exp_times_startSEXP, SEXP indiv_poss_exp_times_endSEXP, SEXP timevarying_groupsSEXP, SEXP variant_specific_parsSEXP, SEXP tempSEXP, SEXP solve_likelihoodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type unique_theta_indices(unique_theta_indicesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type unique_biomarker_groups(unique_biomarker_groupsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type indiv_group_indices(indiv_group_indicesSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type infection_history_mat(infection_history_matSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type infection_history_mat_indices(infection_history_mat_indicesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type likelihoods_pre_proposal(likelihoods_pre_proposalSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type sampled_indivs(sampled_indivsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type n_times_samp_vec(n_times_samp_vecSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type age_mask(age_maskSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type sample_mask(sample_maskSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type n_alive(n_aliveSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type n_infections(n_infectionsSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type n_infected_group(n_infected_groupSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type prior_lookup(prior_lookupSEXP);
    Rcpp::traits::input_parameter< const double& >::type proposal_inf_hist_indiv_swap_ratio(proposal_inf_hist_indiv_swap_ratioSEXP);
    Rcpp::traits::input_parameter< const int& >::type swap_distance(swap_distanceSEXP);
    Rcpp::traits::input_parameter< const bool& >::type propose_from_prior(propose_from_priorSEXP);
    Rcpp::traits::input_parameter< const double& >::type shape1(shape1SEXP);
    Rcpp::traits::input_parameter< const double& >::type shape2(shape2SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type possible_exposure_times(possible_exposure_timesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type possible_exposure_times_indices(possible_exposure_times_indicesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type exposure_groups(exposure_groupsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type unique_exposure_groups(unique_exposure_groupsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sample_times(sample_timesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type type_data_start(type_data_startSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type biomarker_groups(biomarker_groupsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type sample_data_start(sample_data_startSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type antibody_data_start(antibody_data_startSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type nrows_per_sample(nrows_per_sampleSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type cum_nrows_per_individual_in_data(cum_nrows_per_individual_in_dataSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type cum_nrows_per_individual_in_repeat_data(cum_nrows_per_individual_in_repeat_dataSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type popn_group_id_vec(popn_group_id_vecSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type biomarker_id_indices(biomarker_id_indicesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type start_level_indices(start_level_indicesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type starting_antibody_levels(starting_antibody_levelsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type births(birthsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type antigenic_map_long(antigenic_map_longSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type antigenic_map_short(antigenic_map_shortSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type antigenic_distances(antigenic_distancesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type antibody_data(antibody_dataSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type antibody_data_repeats(antibody_data_repeatsSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_measurements_total(n_measurements_totalSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type repeat_indices(repeat_indicesSEXP);
    Rcpp::traits::input_parameter< const bool& >::type repeat_data_exist(repeat_data_existSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type measurement_shifts(measurement_shiftsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type proposal_iter(proposal_iterSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type accepted_iter(accepted_iterSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type proposal_swap(proposal_swapSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type accepted_swap(accepted_swapSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type overall_swap_proposals(overall_swap_proposalsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type overall_add_proposals(overall_add_proposalsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type time_sample_probs(time_sample_probsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type total_alive(total_aliveSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type data_types(data_typesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type obs_weights(obs_weightsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type indiv_possible_exposure_times_indices(indiv_possible_exposure_times_indicesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type indiv_poss_exp_times_start(indiv_poss_exp_times_startSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type indiv_poss_exp_times_end(indiv_poss_exp_times_endSEXP);
    Rcpp::traits::input_parameter< const bool >::type timevarying_groups(timevarying_groupsSEXP);
    Rcpp::traits::input_parameter< const bool >::type variant_specific_pars(variant_specific_parsSEXP);
    Rcpp::traits::input_parameter< const double >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< bool >::type solve_likelihood(solve_likelihoodSEXP);
    rcpp_result_gen = Rcpp::wrap(inf_hist_prop_prior_v2_and_v4(theta, unique_theta_indices, unique_biomarker_groups, indiv_group_indices, infection_history_mat, infection_history_mat_indices, likelihoods_pre_proposal, sampled_indivs, n_times_samp_vec, age_mask, sample_mask, n_alive, n_infections, n_infected_group, prior_lookup, proposal_inf_hist_indiv_swap_ratio, swap_distance, propose_from_prior, shape1, shape2, possible_exposure_times, possible_exposure_times_indices, exposure_groups, unique_exposure_groups, sample_times, type_data_start, biomarker_groups, sample_data_start, antibody_data_start, nrows_per_sample, cum_nrows_per_individual_in_data, cum_nrows_per_individual_in_repeat_data, popn_group_id_vec, biomarker_id_indices, start_level_indices, starting_antibody_levels, births, antigenic_map_long, antigenic_map_short, antigenic_distances, antibody_data, antibody_data_repeats, n_measurements_total, repeat_indices, repeat_data_exist, measurement_shifts, proposal_iter, accepted_iter, proposal_swap, accepted_swap, overall_swap_proposals, overall_add_proposals, time_sample_probs, total_alive, data_types, obs_weights, indiv_possible_exposure_times_indices, indiv_poss_exp_times_start, indiv_poss_exp_times_end, timevarying_groups, variant_specific_pars, temp, solve_likelihood));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_serosolver_antibody_model", (DL_FUNC) &_serosolver_antibody_model, 26},
    {"_serosolver_subset_nullable_vector", (DL_FUNC) &_serosolver_subset_nullable_vector, 3},
    {"_serosolver_logistic_transform", (DL_FUNC) &_serosolver_logistic_transform, 1},
    {"_serosolver_logit_transform", (DL_FUNC) &_serosolver_logit_transform, 1},
    {"_serosolver_transform_parameters_cpp", (DL_FUNC) &_serosolver_transform_parameters_cpp, 6},
    {"_serosolver_get_starting_antibody_levels", (DL_FUNC) &_serosolver_get_starting_antibody_levels, 3},
    {"_serosolver_sum_likelihoods", (DL_FUNC) &_serosolver_sum_likelihoods, 3},
    {"_serosolver_create_cross_reactivity_vector", (DL_FUNC) &_serosolver_create_cross_reactivity_vector, 2},
    {"_serosolver_sum_buckets", (DL_FUNC) &_serosolver_sum_buckets, 2},
    {"_serosolver_sum_infections_by_group", (DL_FUNC) &_serosolver_sum_infections_by_group, 4},
    {"_serosolver_add_measurement_shifts", (DL_FUNC) &_serosolver_add_measurement_shifts, 4},
    {"_serosolver_inf_mat_prior_cpp", (DL_FUNC) &_serosolver_inf_mat_prior_cpp, 4},
    {"_serosolver_inf_mat_prior_cpp_vector", (DL_FUNC) &_serosolver_inf_mat_prior_cpp_vector, 4},
    {"_serosolver_inf_mat_prior_group_cpp", (DL_FUNC) &_serosolver_inf_mat_prior_group_cpp, 4},
    {"_serosolver_inf_mat_prior_group_cpp_vector", (DL_FUNC) &_serosolver_inf_mat_prior_group_cpp_vector, 4},
    {"_serosolver_inf_mat_prior_total_group_cpp", (DL_FUNC) &_serosolver_inf_mat_prior_total_group_cpp, 4},
    {"_serosolver_likelihood_func_fast", (DL_FUNC) &_serosolver_likelihood_func_fast, 3},
    {"_serosolver_likelihood_func_fast_continuous", (DL_FUNC) &_serosolver_likelihood_func_fast_continuous, 3},
    {"_serosolver_likelihood_func_fast_continuous_fp", (DL_FUNC) &_serosolver_likelihood_func_fast_continuous_fp, 3},
    {"_serosolver_inf_hist_prop_prior_v3", (DL_FUNC) &_serosolver_inf_hist_prop_prior_v3, 10},
    {"_serosolver_inf_hist_prop_prior_v2_and_v4", (DL_FUNC) &_serosolver_inf_hist_prop_prior_v2_and_v4, 63},
    {NULL, NULL, 0}
};

RcppExport void R_init_serosolver(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
