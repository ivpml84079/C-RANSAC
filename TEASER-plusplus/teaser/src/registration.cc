/**
 * Copyright 2020, Massachusetts Institute of Technology,
 * Cambridge, MA 02139
 * All Rights Reserved
 * Authors: Jingnan Shi, et al. (see THANKS for the full author list)
 * See LICENSE for the license information
 */

#include "teaser/registration.h"

#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <iterator>

//C-RANSAC Modify: CHANGE NEW INCLUDE
#include <time.h> 
#include <fstream>
#include <streambuf>
#include <map>
#include <unistd.h> 

#include "teaser/utils.h"
#include "teaser/graph.h"
#include "teaser/macros.h"

#define NOISE_BOUND 0.05
double PrNoise = 2 * NOISE_BOUND;

int unknownScale = 1;
double STswitch = 0;
int first_time = 1;
double scale_noise = 0;
double translation_noise = 0;
double scale_last_best = 0;
Eigen::Matrix3d rotation_last_best;
Eigen::Vector3d translation_last_best;
double rotation_similar = 0.01;
int local_max_iter = 10;
bool longholi = 0; //stop process when we used all points

//C-RANSAC Modify
void teaser::ScalarTLSEstimator::estimate(const Eigen::RowVectorXd& X,
                                          const Eigen::RowVectorXd& ranges, double* estimate,
                                          Eigen::Matrix<bool, 1, Eigen::Dynamic>* inliers) {
  // check input parameters
  bool dimension_inconsistent = (X.rows() != ranges.rows()) || (X.cols() != ranges.cols());
  if (inliers) {
    dimension_inconsistent |= ((inliers->rows() != 1) || (inliers->cols() != ranges.cols()));
  }
  bool only_one_element = (X.rows() == 1) && (X.cols() == 1);
  assert(!dimension_inconsistent);
  assert(!only_one_element); // TODO: admit a trivial solution

  //*****************C-RANSAC Modify: S T CHANGE FROM HERE**********************
  int N = X.cols();
  if(STswitch == 0){ //solve scale
    //RANSAC confidence level
    STswitch++;
    int last_not_best = 0;
    int best_inliers_count = 0;
    srand( time(NULL) );
    double confidence = 0;
    int iteration = 0;
    if(!first_time)
    {
      iteration++;
      for (size_t j= 0 ;j < N ;++j){
        if(abs(X(j) - scale_last_best) <= ranges(j))
          best_inliers_count++;
      }
      *estimate = scale_last_best;
      confidence = 1.0 - pow(1.0 - ((double)best_inliers_count / (double)N), iteration);
      // std::cout << "!first_time scale confidence: " << confidence << "  iteration: " << iteration << std::endl;
      // std::cout << "scale_last_best: " << scale_last_best << std::endl;
    }
    while(confidence < 0.99)
    {
      iteration++;
      int ran = rand() % N;
      int curr_count = 0;
      for (size_t j= 0 ;j < N ;++j){
        if(abs(X(j) - X(ran)) <= ranges(j))
          curr_count++;
      }
      if(curr_count > best_inliers_count)
      {
        best_inliers_count = curr_count;
        *estimate = X(ran);
        if(!first_time) last_not_best = 1;
      }
      confidence = 1.0 - pow(1.0 - ((double)best_inliers_count / (double)N), iteration);
      // std::cout << "while scale confidence: " << confidence << "  iteration: " << iteration << std::endl;
      // std::cout << "scale_last_best: " << scale_last_best << " X(ran): " << X(ran) << " estimate: " << *estimate << std::endl;
    }
    if (inliers && estimate) {
      *inliers = (X.array() - *estimate).array().abs() <= ranges.array();
      double sum_left = 0, sum_right = 0;
      for (size_t i= 0 ;i < inliers->size() ;++i){
        if(inliers->operator()(i) == true)
        {
          sum_left += 1.0 / (ranges(i) * ranges(i));
          sum_right += X(i) / (ranges(i) * ranges(i));
        }
      }
      if(!std::isnan(sum_right) && !std::isnan(sum_left))
        *estimate = sum_right / sum_left;
    }
    // std::cout << "Scale Iteration Times: " << iteration << std::endl;
  }
  else{  //solve translation
    //The max-stabbing problem
    std::vector<std::pair<double, int>> h;
    if(!first_time)
      h.resize(2*(N+1), std::make_pair(0, 0));
    else
      h.resize(2*N, std::make_pair(0, 0));
    int L = 0;
    for (size_t i= 0 ;i < N ; i++){
      h[L].first = X(i) - translation_noise;
      h[L++].second = i;
      h[L].first = X(i) + translation_noise;
      h[L++].second = i;
    }
    double transAxis = 0;
    if(!first_time)
    {
      if(STswitch == 1) //x
      {
        transAxis = translation_last_best(0);
        STswitch++;
      }
      else if(STswitch == 2) //y
      {
        transAxis = translation_last_best(1);
        STswitch++;
      }
      else  //z
      {
        transAxis = translation_last_best(2);
        if(unknownScale)
          STswitch = 0; //scale
        else
          STswitch = 1; //x
      }
      h[L].first = transAxis - translation_noise;
      h[L++].second = N;
      h[L].first = transAxis + translation_noise;
      h[L++].second = N;
      N++;
    }
    std::sort(h.begin(), h.end(), [](std::pair<double, int> a, std::pair<double, int> b) { return a.first < b.first; });

    std::vector<int> record(N, 0);
    int currLine = 0, bestLine = 0, remainingLine = N;
    double sum_left = 0, sum_right = 0;
    for(int i = 0; i < h.size(); i++)
    {
      double x = X(h[i].second);
      if(!first_time && h[i].second == N - 1)
        x = transAxis;
      if(record[h[i].second] == 0) //new line segment
      {
        sum_left += 1.0 / (translation_noise * translation_noise);
        sum_right += x / (translation_noise * translation_noise);
        currLine++; remainingLine--;
        record[h[i].second] = 1;
      }
      else  //remove line segment
      {
        if(currLine > bestLine)
        {
          bestLine = currLine;
          if(!std::isnan(sum_right) && !std::isnan(sum_left))
            *estimate = sum_right / sum_left;
          else
            *estimate = x;
        }
        sum_left -= 1.0 / (translation_noise * translation_noise);
        sum_right -= x / (translation_noise * translation_noise);
        currLine--;
        record[h[i].second] = 0;
        if(currLine + remainingLine <= bestLine) //early stop
          break;
      }
    }
    if (inliers) {
      for (size_t i = 0; i < X.cols(); ++i) {
        if(abs(X(i) - *estimate) <= translation_noise)
          (*inliers)(0, i) = 1;
      }
    }
  }
}

void teaser::ScalarTLSEstimator::estimate_tiled(const Eigen::RowVectorXd& X,
                                                const Eigen::RowVectorXd& ranges, const int& s,
                                                double* estimate,
                                                Eigen::Matrix<bool, 1, Eigen::Dynamic>* inliers) {
  // check input parameters
  bool dimension_inconsistent = (X.rows() != ranges.rows()) || (X.cols() != ranges.cols());
  if (inliers) {
    dimension_inconsistent |= ((inliers->rows() != 1) || (inliers->cols() != ranges.cols()));
  }
  bool only_one_element = (X.rows() == 1) && (X.cols() == 1);
  assert(!dimension_inconsistent);
  assert(!only_one_element); // TODO: admit a trivial solution

  // Prepare variables for calculations
  int N = X.cols();
  Eigen::RowVectorXd h(N * 2);
  h << X - ranges, X + ranges;
  // ascending order
  std::sort(h.data(), h.data() + h.cols(), [](double a, double b) { return a < b; });
  // calculate interval centers
  Eigen::RowVectorXd h_centers = (h.head(h.cols() - 1) + h.tail(h.cols() - 1)) / 2;
  auto nr_centers = h_centers.cols();

  // calculate weights
  Eigen::RowVectorXd weights = ranges.array().square();
  weights = weights.array().inverse();

  Eigen::RowVectorXd x_hat = Eigen::MatrixXd::Zero(1, nr_centers);
  Eigen::RowVectorXd x_cost = Eigen::MatrixXd::Zero(1, nr_centers);

  // loop tiling
  size_t ih_bound = ((nr_centers) & ~((s)-1));
  size_t jh_bound = ((N) & ~((s)-1));

  std::vector<double> ranges_inverse_sum_vec(nr_centers, 0);
  std::vector<double> dot_X_weights_vec(nr_centers, 0);
  std::vector<double> dot_weights_consensus_vec(nr_centers, 0);
  std::vector<std::vector<double>> X_consensus_table(nr_centers, std::vector<double>());

  auto inner_loop_f = [&](const size_t& i, const size_t& jh, const size_t& jl_lower_bound,
                          const size_t& jl_upper_bound) {
    double& ranges_inverse_sum = ranges_inverse_sum_vec[i];
    double& dot_X_weights = dot_X_weights_vec[i];
    double& dot_weights_consensus = dot_weights_consensus_vec[i];
    std::vector<double>& X_consensus_vec = X_consensus_table[i];

    size_t j = 0;
    for (size_t jl = jl_lower_bound; jl < jl_upper_bound; ++jl) {
      j = jh + jl;
      bool consensus = std::abs(X(j) - h_centers(i)) <= ranges(j);
      if (consensus) {
        dot_X_weights += X(j) * weights(j);
        dot_weights_consensus += weights(j);
        X_consensus_vec.push_back(X(j));
      } else {
        ranges_inverse_sum += ranges(j);
      }
    }

    if (j == N - 1) {
      // x_hat(i) = dot(X(consensus), weights(consensus)) / dot(weights, consensus);
      x_hat(i) = dot_X_weights / dot_weights_consensus;

      // residual = X(consensus)-x_hat(i);
      Eigen::Map<Eigen::VectorXd> X_consensus(X_consensus_vec.data(), X_consensus_vec.size());
      Eigen::VectorXd residual = X_consensus.array() - x_hat(i);

      // x_cost(i) = dot(residual,residual) + sum(ranges(~consensus));
      x_cost(i) = residual.squaredNorm() + ranges_inverse_sum;
    }
  };

#pragma omp parallel for default(none) shared(                                                     \
    jh_bound, ih_bound, ranges_inverse_sum_vec, dot_X_weights_vec, dot_weights_consensus_vec,      \
    X_consensus_table, h_centers, weights, N, X, x_hat, x_cost, s, ranges, inner_loop_f)
  for (size_t ih = 0; ih < ih_bound; ih += s) {
    for (size_t jh = 0; jh < jh_bound; jh += s) {
      for (size_t il = 0; il < s; ++il) {
        size_t i = ih + il;
        inner_loop_f(i, jh, 0, s);
      }
    }
  }

  // finish the left over entries
  // 1. Finish the unfinished js
#pragma omp parallel for default(none)                                                             \
    shared(jh_bound, ih_bound, ranges_inverse_sum_vec, dot_X_weights_vec,                          \
           dot_weights_consensus_vec, X_consensus_table, h_centers, weights, N, X, x_hat, x_cost,  \
           s, ranges, nr_centers, inner_loop_f)
  for (size_t i = 0; i < nr_centers; ++i) {
    inner_loop_f(i, 0, jh_bound, N);
  }

  // 2. Finish the unfinished is
#pragma omp parallel for default(none)                                                             \
    shared(jh_bound, ih_bound, ranges_inverse_sum_vec, dot_X_weights_vec,                          \
           dot_weights_consensus_vec, X_consensus_table, h_centers, weights, N, X, x_hat, x_cost,  \
           s, ranges, nr_centers, inner_loop_f)
  for (size_t i = ih_bound; i < nr_centers; ++i) {
    inner_loop_f(i, 0, 0, N);
  }

  size_t min_idx;
  x_cost.minCoeff(&min_idx);
  double estimate_temp = x_hat(min_idx);
  if (estimate) {
    // update estimate output if it's not nullptr
    *estimate = estimate_temp;
  }
  if (inliers) {
    // update inlier output if it's not nullptr
    *inliers = (X.array() - estimate_temp).array().abs() <= ranges.array();
  }
}

void teaser::FastGlobalRegistrationSolver::solveForRotation(
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& src,
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& dst, Eigen::Matrix3d* rotation,
    Eigen::Matrix<bool, 1, Eigen::Dynamic>* inliers) {
  assert(rotation);                 // make sure R is not a nullptr
  assert(src.cols() == dst.cols()); // check dimensions of input data
  assert(params_.gnc_factor > 1);   // make sure mu will decrease
  assert(params_.noise_bound != 0); // make sure noise bound is not zero
  if (inliers) {
    assert(inliers->cols() == src.cols());
  }

  // Prepare some intermediate variables
  double noise_bound_sq = std::pow(params_.noise_bound, 2);
  size_t match_size = src.cols();
  cost_ = std::numeric_limits<double>::infinity();

  // Calculate the initial mu
  double src_diameter = teaser::utils::calculateDiameter<double, 3>(src);
  double dest_diameter = teaser::utils::calculateDiameter<double, 3>(dst);
  double global_scale = src_diameter > dest_diameter ? src_diameter : dest_diameter;
  global_scale /= noise_bound_sq;
  double mu = std::pow(global_scale, 2) / noise_bound_sq;

  // stopping condition for mu
  double min_mu = 1.0;
  *rotation = Eigen::Matrix3d::Identity(3, 3); // rotation matrix
  Eigen::Matrix<double, 1, Eigen::Dynamic> l_pq(1, match_size);
  l_pq.setOnes(1, match_size);

  // Assumptions of the two inputs:
  // they should be of the same scale,
  // outliers should be removed as much as possible
  // input vectors should contain TIM vectors (if only estimating rotation)
  for (size_t i = 0; i < params_.max_iterations; ++i) {
    double scaled_mu = mu * noise_bound_sq;

    // 1. Optimize for line processes weights
    Eigen::Matrix<double, 3, 1> q, p, rpq;
    for (size_t j = 0; j < match_size; ++j) {
      // p = Rq
      q = src.col(j);
      p = dst.col(j);
      rpq = p - (*rotation) * q;
      l_pq(j) = std::pow(scaled_mu / (scaled_mu + rpq.squaredNorm()), 2);
    }

    // 2. Optimize for Rotation Matrix
    *rotation = teaser::utils::svdRot(src, dst, l_pq);

    // update cost
    Eigen::Matrix<double, 3, Eigen::Dynamic> diff = (dst - (*rotation) * src).array().square();
    cost_ = ((scaled_mu * diff.colwise().sum()).array() /
             (scaled_mu + diff.colwise().sum().array()).array())
                .sum();

    // additional termination conditions
    if (cost_ < params_.cost_threshold || mu < min_mu) {
      // TEASER_DEBUG_INFO_MSG("Convergence condition met.");
      // TEASER_DEBUG_INFO_MSG("Iterations: " << i);
      // TEASER_DEBUG_INFO_MSG("Mu: " << mu);
      // TEASER_DEBUG_INFO_MSG("Cost: " << cost_);
      break;
    }

    // update mu
    mu /= params_.gnc_factor;
  }

  if (inliers) {
    *inliers = l_pq.cast<bool>();
  }
}

//params.estimate_scaling = true;
void teaser::TLSScaleSolver::solveForScale(const Eigen::Matrix<double, 3, Eigen::Dynamic>& src,
                                           const Eigen::Matrix<double, 3, Eigen::Dynamic>& dst,
                                           double* scale,
                                           Eigen::Matrix<bool, 1, Eigen::Dynamic>* inliers) {

  Eigen::Matrix<double, 1, Eigen::Dynamic> v1_dist =
      src.array().square().colwise().sum().array().sqrt();
  Eigen::Matrix<double, 1, Eigen::Dynamic> v2_dist =
      dst.array().square().colwise().sum().array().sqrt();

  Eigen::Matrix<double, 1, Eigen::Dynamic> raw_scales = v2_dist.array() / v1_dist.array();
  double beta = 2 * noise_bound_ * sqrt(cbar2_);

  //C-RANSAC Modified
  scale_noise = 2 * noise_bound_ * sqrt(cbar2_);
  Eigen::Matrix<double, 1, Eigen::Dynamic> alphas = beta * v1_dist.cwiseInverse();
  
  tls_estimator_.estimate(raw_scales, alphas, scale, inliers);
}

//params.estimate_scaling = false;
void teaser::ScaleInliersSelector::solveForScale(
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& src,
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& dst, double* scale,
    Eigen::Matrix<bool, 1, Eigen::Dynamic>* inliers) {
  // We assume no scale difference between the two vectors of points.
  *scale = 1;

  Eigen::Matrix<double, 1, Eigen::Dynamic> v1_dist =
      src.array().square().colwise().sum().array().sqrt();
  Eigen::Matrix<double, 1, Eigen::Dynamic> v2_dist =
      dst.array().square().colwise().sum().array().sqrt();
  double beta = 2 * noise_bound_ * sqrt(cbar2_);

  // A pair-wise correspondence is an inlier if it passes the following test:
  // abs(|dst| - |src|) is within maximum allowed error
  *inliers = (v1_dist.array() - v2_dist.array()).array().abs() <= beta;
}

void teaser::TLSTranslationSolver::solveForTranslation(
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& src,
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& dst, Eigen::Vector3d* translation,
    Eigen::Matrix<bool, 1, Eigen::Dynamic>* inliers) {
  assert(src.cols() == dst.cols());
  if (inliers) {
    assert(inliers->cols() == src.cols());
  }

  // Raw translation
  Eigen::Matrix<double, 3, Eigen::Dynamic> raw_translation = dst - src;

  // Error bounds for each measurements
  int N = src.cols();
  double beta = noise_bound_ * sqrt(cbar2_);
  translation_noise = noise_bound_ * sqrt(cbar2_);
  Eigen::Matrix<double, 1, Eigen::Dynamic> alphas = beta * Eigen::MatrixXd::Ones(1, N);

  // Estimate x, y, and z component of translation: perform TLS on each row
  *inliers = Eigen::Matrix<bool, 1, Eigen::Dynamic>::Ones(1, N);
  Eigen::Matrix<bool, 1, Eigen::Dynamic> inliers_temp(1, N);
  for (size_t i = 0; i < raw_translation.rows(); ++i) {
    tls_estimator_.estimate(raw_translation.row(i), alphas, &((*translation)(i)), &inliers_temp);
    // element-wise AND using component-wise product (Eigen 3.2 compatible)
    // a point is an inlier iff. x,y,z are all inliers
    *inliers = (*inliers).cwiseProduct(inliers_temp);
  }
}

teaser::RobustRegistrationSolver::RobustRegistrationSolver(
    const teaser::RobustRegistrationSolver::Params& params) {
  reset(params);
}

Eigen::Matrix<double, 3, Eigen::Dynamic>
teaser::RobustRegistrationSolver::computeTIMs(const Eigen::Matrix<double, 3, Eigen::Dynamic>& v,
                                              Eigen::Matrix<int, 2, Eigen::Dynamic>* map) {

  auto N = v.cols();
  Eigen::Matrix<double, 3, Eigen::Dynamic> vtilde(3, N * (N - 1) / 2);
  map->resize(2, N * (N - 1) / 2);

#pragma omp parallel for default(none) shared(N, v, vtilde, map)
  for (size_t i = 0; i < N - 1; i++) {
    // Calculate some important indices
    // For each measurement, we compute the TIMs between itself and all the measurements after it.
    // For example:
    // i=0: add N-1 TIMs
    // i=1: add N-2 TIMs
    // etc..
    // i=k: add N-1-k TIMs
    // And by arithmatic series, we can get the starting index of each segment be:
    // k*N - k*(k+1)/2
    size_t segment_start_idx = i * N - i * (i + 1) / 2;
    size_t segment_cols = N - 1 - i;

    // calculate TIM
    Eigen::Matrix<double, 3, 1> m = v.col(i);
    Eigen::Matrix<double, 3, Eigen::Dynamic> temp = v - m * Eigen::MatrixXd::Ones(1, N);

    // concatenate to the end of the tilde vector
    vtilde.middleCols(segment_start_idx, segment_cols) = temp.rightCols(segment_cols);

    // populate the index map
    Eigen::Matrix<int, 2, Eigen::Dynamic> map_addition(2, N);
    for (size_t j = 0; j < N; ++j) {
      map_addition(0, j) = i;
      map_addition(1, j) = j;
    }
    map->middleCols(segment_start_idx, segment_cols) = map_addition.rightCols(segment_cols);
  }

  return vtilde;
}

teaser::RegistrationSolution
teaser::RobustRegistrationSolver::solve(const teaser::PointCloud& src_cloud,
                                        const teaser::PointCloud& dst_cloud,
                                        const std::vector<std::pair<int, int>> correspondences) {
  Eigen::Matrix<double, 3, Eigen::Dynamic> src(3, correspondences.size());
  Eigen::Matrix<double, 3, Eigen::Dynamic> dst(3, correspondences.size());
  for (size_t i = 0; i < correspondences.size(); ++i) {
    auto src_idx = std::get<0>(correspondences[i]);
    auto dst_idx = std::get<1>(correspondences[i]);
    src.col(i) << src_cloud[src_idx].x, src_cloud[src_idx].y, src_cloud[src_idx].z;
    dst.col(i) << dst_cloud[dst_idx].x, dst_cloud[dst_idx].y, dst_cloud[dst_idx].z;
  }
  return solve(src, dst);
}

//C-RANSAC Modify
teaser::RegistrationSolution
teaser::RobustRegistrationSolver::solve(const Eigen::Matrix<double, 3, Eigen::Dynamic>& src,
                                        const Eigen::Matrix<double, 3, Eigen::Dynamic>& dst) {
  assert(scale_solver_ && rotation_solver_ && translation_solver_);

  // Handle deprecated params
  if (!params_.use_max_clique) {
    TEASER_DEBUG_INFO_MSG(
        "Using deprecated param field use_max_clique. Switch to inlier_selection_mode instead.");
    params_.inlier_selection_mode = INLIER_SELECTION_MODE::NONE;
  }
  if (!params_.max_clique_exact_solution) {
    TEASER_DEBUG_INFO_MSG("Using deprecated param field max_clique_exact_solution. Switch to "
                          "inlier_selection_mode instead.");
    params_.inlier_selection_mode = INLIER_SELECTION_MODE::PMC_HEU;
  }

  /**
   * Steps to estimate T/R/s
   *
   * Estimate Scale
   * -- compute TIMs
   *
   * Remove outliers
   * De-scale the TIMs
   *        v2tilde = v2tilde/s_est; % correct scale from v2 side, more stable
   *
   * Estimate rotation
   *
   * Estimate Translation
   */

  //*****************C-RANSAC Modify: CHANGE FROM HERE**********************    
  unknownScale = params_.estimate_scaling;
  double scale_best_sampled = 1.0, scale_best_host = 1.0;
  Eigen::Matrix3d rotation_best_sampled, rotation_best_host;
  Eigen::Vector3d translation_best_sampled, translation_best_host;
  rotation_best_sampled << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  rotation_best_host << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  translation_best_sampled << 0, 0, 0;
  translation_best_host << 0, 0, 0;
  scale_last_best = 1.0;
  rotation_last_best << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  translation_last_best << 0, 0, 0;
  int C = src.cols();
  // std::cout << std::endl << std::endl << "Point Set: " << C << std::endl;   
  //L set & L reduced set construct
  Eigen::Matrix<double, 3, Eigen::Dynamic> src_tims_ttgg(3, C * (C - 1) / 2);
  Eigen::Matrix<double, 3, Eigen::Dynamic> dst_tims_ttgg(3, C * (C - 1) / 2);
  Eigen::Matrix<int, 2, Eigen::Dynamic> src_tims_map_ttgg(2, C * (C - 1) / 2);
  Eigen::Matrix<int, 2, Eigen::Dynamic> dst_tims_map_ttgg(2, C * (C - 1) / 2);
  Eigen::Matrix<double, 1, Eigen::Dynamic> X(1, C * (C - 1) / 2);
  int MaxScale = 10000;
  int binsize = 20; //20: 0.05; 40: 0.025;
  std::vector<std::vector<int>> H(MaxScale*binsize, std::vector<int>());
  int max_H_index = 0, max_H_height = 0;
  int L = 0;
  for(int i = 0; i < C - 1; i++)
  {
    for(int j = i + 1; j < C; j++)
    {
      Eigen::Vector3d tempSrc = src.col(j) - src.col(i);
      Eigen::Vector3d tempDst = dst.col(j) - dst.col(i);

      src_tims_ttgg.col(L) = tempSrc;
      dst_tims_ttgg.col(L) = tempDst;
      src_tims_map_ttgg(0, L) = i;
      src_tims_map_ttgg(1, L) = j;
      dst_tims_map_ttgg(0, L) = i;
      dst_tims_map_ttgg(1, L) = j;
      X(L) = tempDst.norm() / tempSrc.norm();
      //Histogram
      if(X(L) > MaxScale)
      {
        MaxScale = ceil(MaxScale + X(L));
        H.resize(MaxScale * binsize, std::vector<int>());
      }
      int H_index = floor((X(L) - 0) / MaxScale * H.size());
      if(H_index == H.size())
        H_index--; //incase of outrange
      else if(H_index > H.size() || H_index < 0)
        H_index = 0;
      H[H_index].push_back(L);
      if(H[H_index].size() > max_H_height)
      {
        max_H_height = H[H_index].size();
        max_H_index = H_index;
      }
      L++;
    }
  }
  Eigen::Matrix<double, 3, Eigen::Dynamic> src_tims_tt(3, L);
  Eigen::Matrix<double, 3, Eigen::Dynamic> dst_tims_tt(3, L);
  Eigen::Matrix<int, 2, Eigen::Dynamic> src_tims_map_tt(2, L);
  Eigen::Matrix<int, 2, Eigen::Dynamic> dst_tims_map_tt(2, L);
  src_tims_tt = src_tims_ttgg.block(0,0,3,L); 
  dst_tims_tt = dst_tims_ttgg.block(0,0,3,L); 
  src_tims_map_tt = src_tims_map_ttgg.block(0,0,3,L); 
  dst_tims_map_tt = dst_tims_map_ttgg.block(0,0,3,L); 
  // std::cout << "Line Vector Set: " << L << std::endl; 
  std::vector<int> L_reduced_set;
  if(params_.estimate_scaling)
  {
    L_reduced_set = H[max_H_index];
    if(max_H_index != 0)
      L_reduced_set.insert(L_reduced_set.end(), H[max_H_index - 1].begin(), H[max_H_index - 1].end());
    if(max_H_index != H.size() - 1)
      L_reduced_set.insert(L_reduced_set.end(), H[max_H_index + 1].begin(), H[max_H_index + 1].end());
    STswitch = 0;
  }
  else
  {
    solveForScale(src_tims_tt, dst_tims_tt); 
    //line inlier after scale
    int scaleInlier = 0;
    for (size_t i = 0; i < scale_inliers_mask_.cols(); ++i)
      if (scale_inliers_mask_(0, i))
        scaleInlier++;
    L_reduced_set.resize(scaleInlier);
    int LL = 0;
    for (size_t i = 0; i < scale_inliers_mask_.cols(); ++i)
      if (scale_inliers_mask_(0, i))
        L_reduced_set[LL++] = i;
    STswitch = 1;
  }
  int N = X.cols();
  // std::cout << "L reduced Set: "<< L_reduced_set.size() << std::endl;
  //L reduced set complete
  int best_inliers_count_host = 0, host_r = 0;
  double pro_host = 0.0, Tpro_host = 0.99;
  bool pro_host_not_over = true;
  srand( time(NULL) );   
  int localCount = 1, basicCount = 1;
  double L_sampled_rate = 0.1; //sampled rate from reduced set
  double b_sampled_rate = 0.3; //basic rate from sampled set
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  while(pro_host_not_over) //RANSAC host start, produced RANSAC local
  {
    //L sampled set construct
    int sampled_first_time = 1;
    int L_sampled_set_size = floor(L_reduced_set.size() * L_sampled_rate); //{L_sampled_rate}% of reduced set
    int *L_sampled_set; //save for X index
    if(L_sampled_set_size == 0)
    {
      L_sampled_set_size = L_reduced_set.size();
      L_sampled_set = new int [L_sampled_set_size];
      for(int i = 0; i < L_sampled_set_size; i++)
      {
        L_sampled_set[i] = L_reduced_set[i];
      }
    }
    else
    {
      L_sampled_set = new int [L_sampled_set_size];
      int *duplicate = new int [L_reduced_set.size()](); //check duplicate
      for(int i = 0; i < L_sampled_set_size; i++)
      {
        int ran = 0;
        do
        {
          ran = rand() % L_reduced_set.size();
        } while (duplicate[ran] == 1);
        duplicate[ran] = 1;
        L_sampled_set[i] = L_reduced_set[ran];
      }
      delete [] duplicate;
    }
    Eigen::Matrix<double, 3, Eigen::Dynamic> src_sampled;
    Eigen::Matrix<double, 3, Eigen::Dynamic> dst_sampled;
    Eigen::Matrix<double, 3, Eigen::Dynamic> src_tims_ttt(3, L_sampled_set_size);
    Eigen::Matrix<double, 3, Eigen::Dynamic> dst_tims_ttt(3, L_sampled_set_size);
    Eigen::Matrix<int, 2, Eigen::Dynamic> src_tims_map_ttt(2, L_sampled_set_size);
    Eigen::Matrix<int, 2, Eigen::Dynamic> dst_tims_map_ttt(2, L_sampled_set_size);
    int* dub = new int[src.cols()]();
    for(int i = 0; i < L_sampled_set_size; i++)
    {
      src_tims_ttt.col(i) = src_tims_tt.col(L_sampled_set[i]);
      dst_tims_ttt.col(i) = dst_tims_tt.col(L_sampled_set[i]);
      src_tims_map_ttt.col(i) = src_tims_map_tt.col(L_sampled_set[i]);
      dst_tims_map_ttt.col(i) = dst_tims_map_tt.col(L_sampled_set[i]);
      if(dub[src_tims_map_tt(0, L_sampled_set[i])] == 0)
      {
        src_sampled.conservativeResize(src_sampled.rows(), src_sampled.cols()+1);
        src_sampled.col(src_sampled.cols()-1) = src.col(src_tims_map_tt(0, L_sampled_set[i]));
        dst_sampled.conservativeResize(dst_sampled.rows(), dst_sampled.cols()+1);
        dst_sampled.col(dst_sampled.cols()-1) = dst.col(dst_tims_map_tt(0, L_sampled_set[i]));
        dub[src_tims_map_tt(0, L_sampled_set[i])] = 1;
      }
      if(dub[src_tims_map_tt(1, L_sampled_set[i])] == 0)
      {
        src_sampled.conservativeResize(src_sampled.rows(), src_sampled.cols()+1);
        src_sampled.col(src_sampled.cols()-1) = src.col(src_tims_map_tt(1, L_sampled_set[i]));
        dst_sampled.conservativeResize(dst_sampled.rows(), dst_sampled.cols()+1);
        dst_sampled.col(dst_sampled.cols()-1) = dst.col(dst_tims_map_tt(1, L_sampled_set[i]));
        dub[src_tims_map_tt(1, L_sampled_set[i])] = 1;
      }
    }
    delete [] dub;
    //L sampled set complete

    int best_inliers_count_sampled = 0, local_r = 0;
    double pro_local = 0, Tpro_local = 0.99;
    int pro_local_times = 1;
    bool pro_local_not_over = true;
    int curr_count_zero_count = 0;
    int *duplicate1 = new int [L_sampled_set_size]();
    while (pro_local_not_over) //pro_local calculate start
    {
      // std::cout << "\n*****RANSAC local " << localCount << "  basic " << basicCount << "*****\n";
      // std::cout << "L sampled Set: "<< floor(L_reduced_set.size() * 0.1) << std::endl;
      // std::cout << "src_sampled: " << src_sampled.cols() << std::endl;
      int basic_choose = L_sampled_set_size * b_sampled_rate;
      src_tims_.resize(3, basic_choose);
      dst_tims_.resize(3, basic_choose);
      src_tims_map_.resize(2, basic_choose);
      dst_tims_map_.resize(2, basic_choose);
      int* dub2 = new int [L_sampled_set_size]();
      int x = 0;        
        int dup1 = 0;
      for(int i = 0; i < basic_choose; i++)
      {
        do
        {
          x = rand() % L_sampled_set_size;
        } while (dub2[x]);
        src_tims_.col(i) = src_tims_ttt.col(x);
        dst_tims_.col(i) = dst_tims_ttt.col(x);
        src_tims_map_.col(i) = src_tims_map_ttt.col(x);
        dst_tims_map_.col(i) = dst_tims_map_ttt.col(x);        

        if(duplicate1[x] == 1)
          dup1++;
        duplicate1[x] = 1;

        dub2[x] = 1;
      }
      delete []dub2;

      // std::cout << "basic_choose: " << basic_choose << std::endl;

      reset(params_); //reset        
      params_.noise_bound = 0.05;
      params_.cbar2 = 1;
      params_.estimate_scaling = unknownScale;
      params_.rotation_max_iterations = 100;
      params_.rotation_gnc_factor = 1.4;
      params_.rotation_estimation_algorithm =
          teaser::RobustRegistrationSolver::ROTATION_ESTIMATION_ALGORITHM::GNC_TLS;
      params_.rotation_cost_threshold = 0.005;
      // Handle deprecated params
      if (!params_.use_max_clique) {
        TEASER_DEBUG_INFO_MSG(
            "Using deprecated param field use_max_clique. Switch to inlier_selection_mode instead.");
        params_.inlier_selection_mode = INLIER_SELECTION_MODE::NONE;
      }
      if (!params_.max_clique_exact_solution) {
        TEASER_DEBUG_INFO_MSG("Using deprecated param field max_clique_exact_solution. Switch to "
                              "inlier_selection_mode instead.");
        params_.inlier_selection_mode = INLIER_SELECTION_MODE::PMC_HEU;
      }

      if(params_.estimate_scaling)
      {
        solveForScale(src_tims_, dst_tims_);

        // line inlier after scale
        int scaleInlier = 0;
        for (size_t i = 0; i < scale_inliers_mask_.cols(); ++i)
          if (scale_inliers_mask_(0, i))
            scaleInlier++;
        // std::cout << "Line reduced inliers after scale: " << scaleInlier << std::endl;
        pruned_src_tims_.conservativeResize(pruned_src_tims_.rows(), scaleInlier);
        pruned_dst_tims_.conservativeResize(pruned_dst_tims_.rows(), scaleInlier);
        src_tims_map_rotation_.conservativeResize(src_tims_map_rotation_.rows(), scaleInlier);
        dst_tims_map_rotation_.conservativeResize(dst_tims_map_rotation_.rows(), scaleInlier);
        scaleInlier = 0;
        for (size_t i = 0; i < scale_inliers_mask_.cols(); ++i) {
          if (scale_inliers_mask_(0, i)) {
            pruned_src_tims_.col(scaleInlier) = src_tims_.col(i);
            pruned_dst_tims_.col(scaleInlier) = dst_tims_.col(i);
            src_tims_map_rotation_.col(scaleInlier) = src_tims_map_.col(i);
            dst_tims_map_rotation_.col(scaleInlier++) = dst_tims_map_.col(i);
          }
        }
      }
      else
      {
        solveForScale(src_tims_, dst_tims_);
        pruned_src_tims_ = src_tims_;
        pruned_dst_tims_ = dst_tims_;
        src_tims_map_rotation_ = src_tims_map_;
        dst_tims_map_rotation_ = dst_tims_map_;
      }
      // std::cout << "Line basic size: " << src_tims_.cols() << std::endl;

      // std::cout << std::endl << "host r: " << host_r << " local r: " << local_r << std::endl;
      //maximum clique      
      Eigen::Matrix<double, 3, Eigen::Dynamic> src_inliers_ot;
      Eigen::Matrix<double, 3, Eigen::Dynamic> dst_inliers_ot;
      if(b_sampled_rate == 1.0) // using maximum clique
      {
        if (params_.inlier_selection_mode != INLIER_SELECTION_MODE::NONE) {
          inlier_graph_.populateVertices(src.cols());
          for (size_t i = 0; i < scale_inliers_mask_.cols(); ++i) {
            if (scale_inliers_mask_(0, i)) {
              inlier_graph_.addEdge(src_tims_map_(0, i), src_tims_map_(1, i));
            }
          }

          teaser::MaxCliqueSolver::Params clique_params;

          if (params_.inlier_selection_mode == INLIER_SELECTION_MODE::PMC_EXACT) {
            clique_params.solver_mode = teaser::MaxCliqueSolver::CLIQUE_SOLVER_MODE::PMC_EXACT;
          } else if (params_.inlier_selection_mode == INLIER_SELECTION_MODE::PMC_HEU) {
            clique_params.solver_mode = teaser::MaxCliqueSolver::CLIQUE_SOLVER_MODE::PMC_HEU;
          } else {
            clique_params.solver_mode = teaser::MaxCliqueSolver::CLIQUE_SOLVER_MODE::KCORE_HEU;
          }
          clique_params.time_limit = params_.max_clique_time_limit;
          clique_params.kcore_heuristic_threshold = params_.kcore_heuristic_threshold;

          teaser::MaxCliqueSolver clique_solver(clique_params);
          max_clique_ = clique_solver.findMaxClique(inlier_graph_);
          std::sort(max_clique_.begin(), max_clique_.end());
          TEASER_DEBUG_INFO_MSG("Max Clique of scale estimation inliers: ");
      #ifndef NDEBUG
          std::copy(max_clique_.begin(), max_clique_.end(), std::ostream_iterator<int>(std::cout, " "));
          std::cout << std::endl;
      #endif
          // Abort if max clique size <= 1
          if (max_clique_.size() <= 1) {
            TEASER_DEBUG_INFO_MSG("Clique size too small. Abort.");
            solution_.valid = false;
            return solution_;
          }
        } else {
          // not using clique filtering is equivalent to saying all measurements are in the max clique
          max_clique_.reserve(src.cols());
          for (size_t i = 0; i < src.cols(); ++i) {
            max_clique_.push_back(i);
          }
        }
        // complete graph
        // select the inlier measurements with max clique
        Eigen::Matrix<double, 3, Eigen::Dynamic> src_inliers(3, max_clique_.size());
        Eigen::Matrix<double, 3, Eigen::Dynamic> dst_inliers(3, max_clique_.size());
        std::vector<int> src_inliers_index_to_src(max_clique_.size());
        std::vector<int> dst_inliers_index_to_dst(max_clique_.size());
        for (size_t i = 0; i < max_clique_.size(); ++i) {
          src_inliers.col(i) = src.col(max_clique_[i]);
          dst_inliers.col(i) = dst.col(max_clique_[i]);
          src_inliers_index_to_src[i] = max_clique_[i];
          dst_inliers_index_to_dst[i] = max_clique_[i];
        }
        
        src_inliers_ot = src_inliers;
        dst_inliers_ot = dst_inliers;
        // std::cout << "Line left after maximum clique: " << pruned_src_tims_.cols() << std::endl;
      }

      // Remove scaling for rotation estimation
      pruned_dst_tims_ *= (1 / solution_.scale);
      // Update GNC rotation solver's noise bound with the new information
      // Note: this implicitly assumes that rotation_solver_'s noise bound
      // is set to the original noise bound of the measurements.
      auto params = rotation_solver_->getParams();
      params.noise_bound *= (2 / solution_.scale);
      rotation_solver_->setParams(params);
      // Solve for rotation
      TEASER_DEBUG_INFO_MSG("Starting rotation solver.");
      solveForRotation(pruned_src_tims_, pruned_dst_tims_);
      TEASER_DEBUG_INFO_MSG("Rotation estimation complete.");

      Eigen::Matrix<double, 3, Eigen::Dynamic> rotation_pruned_src(3, 0);
      Eigen::Matrix<double, 3, Eigen::Dynamic> rotation_pruned_dst(3, 0);
      int* dub = new int [src.cols()]();
      for (size_t i = 0; i < rotation_inliers_mask_.cols(); ++i) {
        if (rotation_inliers_mask_[i]) {
          if(src_tims_map_rotation_(0, i) >= src.cols() || src_tims_map_rotation_(0, i) < 0)
            std::cout << "i: " << i << " rotation_inliers_mask_.cols(): " << rotation_inliers_mask_.cols() << " src_tims_map_rotation_(0, i) src: " << src_tims_map_rotation_(0, i) << "  src.cols(): " << src.cols() << std::endl;
          if(dst_tims_map_rotation_(0, i) >= dst.cols() || dst_tims_map_rotation_(0, i) < 0)
            std::cout << "i: " << i << " rotation_inliers_mask_.cols(): " << rotation_inliers_mask_.cols() << " src_tims_map_rotation_(0, i) dst: " << src_tims_map_rotation_(0, i) << "  dst.cols(): " << dst.cols()<< std::endl;
          if(src_tims_map_rotation_(1, i) >= src.cols() || src_tims_map_rotation_(1, i) < 0)
            std::cout << "i: " << i << " rotation_inliers_mask_.cols(): " << rotation_inliers_mask_.cols() << " src_tims_map_rotation_(1, i) src: " << src_tims_map_rotation_(1, i) << "  src.cols(): " << src.cols()<< std::endl;
          if(dst_tims_map_rotation_(1, i) >= dst.cols() || dst_tims_map_rotation_(1, i) < 0)
            std::cout << "i: " << i << " rotation_inliers_mask_.cols(): " << rotation_inliers_mask_.cols() << " src_tims_map_rotation_(1, i) dst: " << src_tims_map_rotation_(1, i) << "  dst.cols(): " << dst.cols()<< std::endl;
          
          if(src_tims_map_rotation_(0, i) >= src.cols() || src_tims_map_rotation_(0, i) < 0)
            continue;
          if(dst_tims_map_rotation_(0, i) >= dst.cols() || dst_tims_map_rotation_(0, i) < 0)
            continue;
          if(src_tims_map_rotation_(1, i) >= src.cols() || src_tims_map_rotation_(1, i) < 0)
            continue;
          if(dst_tims_map_rotation_(1, i) >= dst.cols() || dst_tims_map_rotation_(1, i) < 0)
            continue;
          if(dub[src_tims_map_rotation_(0, i)] == 0)
          {
            rotation_pruned_src.conservativeResize(rotation_pruned_src.rows(), rotation_pruned_src.cols()+1);
            rotation_pruned_src.col(rotation_pruned_src.cols()-1) = src.col(src_tims_map_rotation_(0, i));
            rotation_pruned_dst.conservativeResize(rotation_pruned_dst.rows(), rotation_pruned_dst.cols()+1);
            rotation_pruned_dst.col(rotation_pruned_dst.cols()-1) = dst.col(dst_tims_map_rotation_(0, i));
            dub[src_tims_map_rotation_(0, i)] = 1;
          }
          if(dub[src_tims_map_rotation_(1, i)] == 0)
          {
            rotation_pruned_src.conservativeResize(rotation_pruned_src.rows(), rotation_pruned_src.cols()+1);
            rotation_pruned_src.col(rotation_pruned_src.cols()-1) = src.col(src_tims_map_rotation_(1, i));
            rotation_pruned_dst.conservativeResize(rotation_pruned_dst.rows(), rotation_pruned_dst.cols()+1);
            rotation_pruned_dst.col(rotation_pruned_dst.cols()-1) = dst.col(dst_tims_map_rotation_(1, i));
            dub[src_tims_map_rotation_(1, i)] = 1;
          }
        }
      }
      delete []dub;
      // std::cout << "Point left after rotation: " << rotation_pruned_src.cols() << std::endl;

      if(b_sampled_rate == 1.0) //using points from maximum clique
      {
        rotation_pruned_src = src_inliers_ot;
        rotation_pruned_dst = dst_inliers_ot;
        // std::cout << "Point left after rotation then maximum clique: " << rotation_pruned_src.cols() << std::endl;
      }

      // Solve for translation
      TEASER_DEBUG_INFO_MSG("Starting translation solver.");
      solveForTranslation(solution_.scale * solution_.rotation * rotation_pruned_src, rotation_pruned_dst);
      //Important modified
      solution_.translation /= solution_.scale;
      TEASER_DEBUG_INFO_MSG("Translation estimation complete.");
      // std::cout << "t: \n" << solution_.translation << std::endl;
      if(params_.estimate_scaling)
        STswitch = 0;
      else
        STswitch = 1;
      basicCount++;

      // std::cout << "rotation similar: " << std::abs(std::acos(fmin(fmax(((rotation_last_best.transpose() * solution_.rotation).trace() - 1) / 2, -1.0), 1.0))) << std::endl;
      // std::cout << "translation similar: " << (translation_last_best - solution_.translation).norm() << std::endl;
      if(!first_time && 
          abs(scale_last_best - solution_.scale) <= scale_noise &&
          std::abs(std::acos(fmin(fmax(((rotation_last_best.transpose() * solution_.rotation).trace() - 1) / 2, -1.0), 1.0))) <= rotation_similar &&
          (translation_last_best - solution_.translation).norm() <= translation_noise )
      {
        if(sampled_first_time)
        {
          local_r += host_r + 1;  
          // std::cout << "Solution similar to host.\n";
          // std::cout << "local r = " << local_r << std::endl;
        }
        else
        {
          local_r++;
          // std::cout << "Solution similar to previous basic.\n";
          // std::cout << "local r = " << local_r << std::endl;
        }
        pro_local = 1.0; //local over
        scale_best_sampled = solution_.scale;
        rotation_best_sampled = solution_.rotation;
        translation_best_sampled = solution_.translation;
      }
      else
      {
        local_r++;
        
        if(!first_time && b_sampled_rate < 1.0) //last best inlier rate in RANSAC local(negative effect when using maximumclique)
        {
          // Homogeneous coordinates
          Eigen::Matrix<double, 4, Eigen::Dynamic> src_sampled_h;
          src_sampled_h.resize(4, src_sampled.cols());
          src_sampled_h.topRows(3) = src_sampled;
          src_sampled_h.bottomRows(1) = Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(N);
          Eigen::Matrix4d TRANSFORM;
          TRANSFORM <<  1, 0, 0, 0,
                        0, 1, 0, 0,
                        0, 0, 1, 0,
                        0, 0, 0, 1;
          TRANSFORM.topLeftCorner(3, 3) = rotation_last_best;
          TRANSFORM.topRightCorner(3, 1) = translation_last_best;
          Eigen::Matrix<double, 3, Eigen::Dynamic> src_solve = (scale_last_best * TRANSFORM * src_sampled_h).topRows(3);

          // Find the sampled inliers
          int curr_count = 0;
          for (size_t j= 0 ;j < src_sampled.cols() ;++j){
            if((sqrt(pow(dst_sampled(0, j) - src_solve(0, j), 2) + pow(dst_sampled(1, j) - src_solve(1, j), 2) + pow(dst_sampled(2, j) - src_solve(2, j), 2)) <= PrNoise))
              curr_count++;
          }
          best_inliers_count_sampled = curr_count;
          scale_best_sampled = scale_last_best;
          rotation_best_sampled = rotation_last_best;
          translation_best_sampled = translation_last_best;
        }

        // Homogeneous coordinates
        Eigen::Matrix<double, 4, Eigen::Dynamic> src_sampled_h;
        src_sampled_h.resize(4, src_sampled.cols());
        src_sampled_h.topRows(3) = src_sampled;
        src_sampled_h.bottomRows(1) = Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(N);
        Eigen::Matrix4d TRANSFORM;
        TRANSFORM <<  1, 0, 0, 0,
                      0, 1, 0, 0,
                      0, 0, 1, 0,
                      0, 0, 0, 1;
        TRANSFORM.topLeftCorner(3, 3) = solution_.rotation;
        TRANSFORM.topRightCorner(3, 1) = solution_.translation;
        Eigen::Matrix<double, 3, Eigen::Dynamic> src_solve = (solution_.scale * TRANSFORM * src_sampled_h).topRows(3);

        // Find the sampled inliers
        int curr_count = 0;
        for (size_t j= 0 ;j < src_sampled.cols() ;++j){
          if((sqrt(pow(dst_sampled(0, j) - src_solve(0, j), 2) + pow(dst_sampled(1, j) - src_solve(1, j), 2) + pow(dst_sampled(2, j) - src_solve(2, j), 2)) <= PrNoise))
            curr_count++;
        }
        if(curr_count > best_inliers_count_sampled || first_time)
        {
          // std::cout << "local answer change" << std::endl;
          // std::cout << "curr_count: " << curr_count << std::endl;
          scale_best_sampled = solution_.scale;
          rotation_best_sampled = solution_.rotation;
          translation_best_sampled = solution_.translation;
          best_inliers_count_sampled = curr_count;
        }

        scale_last_best = scale_best_sampled;
        rotation_last_best = rotation_best_sampled;
        translation_last_best = translation_best_sampled;

        pro_local = 1.0 - pow(1.0 - (double)((double)best_inliers_count_sampled / (double)src_sampled.cols()), local_r);
        // std::cout << "local point inliers = " << best_inliers_count_sampled << std::endl;
        // std::cout << "local r = " << local_r << std::endl;
        // std::cout << "Pr local = " << pro_local << std::endl;

        first_time = 0;

        // std::cout << "L_sampled_rate: " << L_sampled_rate << std::endl;
        // std::cout << "b_sampled_rate: " << b_sampled_rate << std::endl;
        if((local_r >= local_max_iter && pro_local <= 0.2) || b_sampled_rate == 1.0) //Need to stop and adjust L_sampled_rate or b_sampled_rate because the current outlier rate is too high
        {
          pro_local = 1.0;
          if(L_sampled_rate == 0.1 && b_sampled_rate == 0.3)
          {
            L_sampled_rate = 0.2; b_sampled_rate = 0.3;
          }
          else if(L_sampled_rate == 0.2 && b_sampled_rate == 0.3)
          {
            L_sampled_rate = 0.5; b_sampled_rate = 0.3;
          }
          else if(L_sampled_rate == 0.5 && b_sampled_rate == 0.3)
          {
            L_sampled_rate = 1.0; b_sampled_rate = 1.0;
          }
        }
      }

      if(pro_local > Tpro_local) //pro_local calculate start
      {
        host_r += local_r;          
        // Homogeneous coordinates
        Eigen::Matrix<double, 4, Eigen::Dynamic> src_h;
        src_h.resize(4, src.cols());
        src_h.topRows(3) = src;
        src_h.bottomRows(1) = Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(N);
        Eigen::Matrix4d TRANSFORM;
        TRANSFORM <<  1, 0, 0, 0,
                      0, 1, 0, 0,
                      0, 0, 1, 0,
                      0, 0, 0, 1;
        TRANSFORM.topLeftCorner(3, 3) = rotation_best_sampled;
        TRANSFORM.topRightCorner(3, 1) = translation_best_sampled;
        Eigen::Matrix<double, 3, Eigen::Dynamic> src_solve = (scale_best_sampled * TRANSFORM * src_h).topRows(3);

        // Find the final inliers
        int curr_count = 0;
        for (size_t j= 0 ;j < src.cols() ;++j){
          if((sqrt(pow(dst(0, j) - src_solve(0, j), 2) + pow(dst(1, j) - src_solve(1, j), 2) + pow(dst(2, j) - src_solve(2, j), 2)) <= PrNoise))
            curr_count++;
        }
        if(curr_count > best_inliers_count_host || pro_host == 0.0 || (b_sampled_rate == 1.0 && curr_count >= best_inliers_count_host))
        {
          // std::cout << "host answer change" << std::endl;
          scale_best_host = scale_best_sampled;
          rotation_best_host = rotation_best_sampled;
          translation_best_host = translation_best_sampled;
          best_inliers_count_host = curr_count;
        }

        scale_last_best = scale_best_host;
        rotation_last_best = rotation_best_host;
        translation_last_best = translation_best_host;

        pro_host = 1.0 - pow(1.0 - (double)((double)best_inliers_count_host / (double)src.cols()), host_r);
        // std::cout << "host point inliers = " << best_inliers_count_host << std::endl;
        // std::cout << "host r = " << host_r << std::endl;
        // std::cout << "Pr host = " << pro_host << std::endl;
        // std::cout << "longholi = " << longholi << std::endl;
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        double curr_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
        if(pro_host > Tpro_host || longholi || curr_time > 60.0)
        {
          pro_host_not_over = false; //stop RANSAC host
          pro_local_not_over = false; //stop RANSAC local
        }
        else
          pro_local_not_over = false; //stop RANSAC local, and go back to RANSAC host

        if(L_sampled_rate == 1.0 && b_sampled_rate == 1.0)
          longholi = true;

        localCount++;
        basicCount = 1;
      }
      sampled_first_time = 0;
    }
    delete [] duplicate1;
    delete [] L_sampled_set;
  }

  solution_.scale = scale_best_host;
  solution_.rotation = rotation_best_host;
  solution_.translation = translation_best_host;

  // Update validity flag
  solution_.valid = true;
  first_time = 1;
  longholi = false;
  return solution_;
}

double teaser::RobustRegistrationSolver::solveForScale(
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& v1,
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& v2) {
  scale_inliers_mask_.resize(1, v1.cols());
  scale_solver_->solveForScale(v1, v2, &(solution_.scale), &scale_inliers_mask_);
  return solution_.scale;
}

Eigen::Vector3d teaser::RobustRegistrationSolver::solveForTranslation(
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& v1,
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& v2) {
  translation_inliers_mask_.resize(1, v1.cols());
  translation_solver_->solveForTranslation(v1, v2, &(solution_.translation),
                                           &translation_inliers_mask_);
  return solution_.translation;
}

Eigen::Matrix3d teaser::RobustRegistrationSolver::solveForRotation(
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& v1,
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& v2) {
  rotation_inliers_mask_.resize(1, v1.cols());
  rotation_solver_->solveForRotation(v1, v2, &(solution_.rotation), &rotation_inliers_mask_);
  return solution_.rotation;
}

//C-RANSAC Modify
void teaser::GNCTLSRotationSolver::solveForRotation(
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& src,
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& dst, Eigen::Matrix3d* rotation,
    Eigen::Matrix<bool, 1, Eigen::Dynamic>* inliers) {
  assert(rotation);                 // make sure R is not a nullptr
  assert(src.cols() == dst.cols()); // check dimensions of input data
  assert(params_.gnc_factor > 1);   // make sure mu will increase
  assert(params_.noise_bound != 0); // make sure noise sigma is not zero
  if (inliers) {
    assert(inliers->cols() == src.cols());
  }

  /**
   * Loop: terminate when:
   *    1. the change in cost in two consecutive runs is smaller than a user-defined threshold
   *    2. # iterations exceeds the maximum allowed
   *
   * Within each loop:
   * 1. fix weights and solve for R
   * 2. fix R and solve for weights
   */

  // Prepare some variables
  size_t match_size = src.cols(); // number of correspondences

  double mu = 1; // arbitrary starting mu

  double prev_cost = std::numeric_limits<double>::infinity();
  cost_ = std::numeric_limits<double>::infinity();
  double noise_bound_sq = std::pow(params_.noise_bound, 2);
  if (noise_bound_sq < 1e-16) {
    noise_bound_sq = 1e-2;
  }
  TEASER_DEBUG_INFO_MSG("GNC rotation estimation noise bound:" << params_.noise_bound);
  TEASER_DEBUG_INFO_MSG("GNC rotation estimation noise bound squared:" << noise_bound_sq);

  Eigen::Matrix<double, 3, Eigen::Dynamic> diffs(3, match_size);
  Eigen::Matrix<double, 1, Eigen::Dynamic> weights(1, match_size);
  Eigen::Matrix<double, 1, Eigen::Dynamic> residuals_sq(1, match_size); 
  
  //************************************C-RANSAC Modify: R change from here****************************************
  int rotation_first_time = 1;  
  //C-RANSAC Modify: Reset between 100 try
  if(first_time == 1)
    rotation_last_best <<  1, 0, 0,
                          0, 1, 0, 
                          0, 0, 1;

  weights.setOnes(1, match_size);
  // Loop for performing GNC-TLS
  for (size_t i = 0; i < params_.max_iterations; ++i) {
    // Fix weights and perform SVD rotation estimation
    // C-RANSAC Modify:
    if(!first_time && rotation_first_time)
    {
      *rotation = rotation_last_best;
      rotation_first_time = 0;
    }
    else
      *rotation = teaser::utils::svdRot(src, dst, weights);

    // Calculate residuals squared
    diffs = (dst - (*rotation) * src).array().square();
    residuals_sq = diffs.colwise().sum();
    if (i == 0) {
      // Initialize rule for mu
      double max_residual = residuals_sq.maxCoeff();
      mu = 1 / (2 * max_residual / noise_bound_sq - 1);
      // Degenerate case: mu = -1 because max_residual is very small
      // i.e., little to none noise
      if (mu <= 0) {
        TEASER_DEBUG_INFO_MSG(
            "GNC-TLS terminated because maximum residual at initialization is very small.");
        break;
      }
    }

    // Fix R and solve for weights in closed form
    double th1 = (mu + 1) / mu * noise_bound_sq;
    double th2 = mu / (mu + 1) * noise_bound_sq;
    cost_ = 0;
    for (size_t j = 0; j < match_size; ++j) {
      // Also calculate cost in this loop
      // Note: the cost calculated is using the previously solved weights
      cost_ += weights(j) * residuals_sq(j);

      if (residuals_sq(j) >= th1) {
        weights(j) = 0;
      } else if (residuals_sq(j) <= th2) {
        weights(j) = 1;
      } else {
        weights(j) = sqrt(noise_bound_sq * mu * (mu + 1) / residuals_sq(j)) - mu;
        assert(weights(j) >= 0 && weights(j) <= 1);
      }
    }

    // Calculate cost
    double cost_diff = std::abs(cost_ - prev_cost);

    // Increase mu
    mu = mu * params_.gnc_factor;
    prev_cost = cost_;

    if (cost_diff < params_.cost_threshold) {
      TEASER_DEBUG_INFO_MSG("GNC-TLS solver terminated due to cost convergence.");
      TEASER_DEBUG_INFO_MSG("Cost diff: " << cost_diff);
      TEASER_DEBUG_INFO_MSG("Iterations: " << i);
      break;
    }
  }

  // C-RANSAC Modify:
  if (inliers) {
    int gf = 0;
    for (size_t i = 0; i < weights.cols(); ++i) {
      if(weights(0, i) >= 0.5)
      {
        (*inliers)(0, i) = 1;
        gf++;
      }
    }
    if(gf <= 10) //Prevent inlier set fail
    {
      for (size_t i = 0; i < weights.cols(); ++i) {
        (*inliers)(0, i) = 1;
      }
    }
  }
}
