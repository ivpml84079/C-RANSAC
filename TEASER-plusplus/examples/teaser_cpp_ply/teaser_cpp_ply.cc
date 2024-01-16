// An example showing TEASER++ registration with the Stanford bunny model
#include <chrono>
#include <iostream>
#include <random>

#include <Eigen/Core>

#include <teaser/ply_io.h>
#include <teaser/registration.h>

//CHANGE NEW INCLUDE 
#include <fstream>
#include <streambuf>
#include <unistd.h> 
#include <cmath>

// Macro constants for generating noise and outliers
#define NOISE_BOUND 0.05
#define N_OUTLIERS_RATE 0.9
#define PI 3.1415926

int unknownScale = 1;

inline double getAngularError(Eigen::Matrix3d R_exp, Eigen::Matrix3d R_est) {
  return std::abs(std::acos(fmin(fmax(((R_exp.transpose() * R_est).trace() - 1) / 2, -1.0), 1.0)))*180/PI;
}

void addNoiseAndOutliers(Eigen::Matrix<double, 3, Eigen::Dynamic>& tgt) {
  Eigen::Matrix<double, 3, Eigen::Dynamic> noise =
      Eigen::Matrix<double, 3, Eigen::Dynamic>::Random(3, tgt.cols()) * NOISE_BOUND;
  NOISE_BOUND / 2;
  tgt = tgt + noise;

  // Add outliers
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis2(0, tgt.cols() - 1); // pos of outliers
  std::uniform_real_distribution<> dis3_lb(-10.0, -5.0); //-10~-5, 5~10
  std::uniform_real_distribution<> dis3_ub(5.0, 10.0);
  std::uniform_real_distribution<> disC(0.0, 1.0);
  std::vector<bool> expected_outlier_mask(tgt.cols(), false);
  int outliers = tgt.cols() * N_OUTLIERS_RATE[N_index];
  for (int i = 0; i < outliers; ++i) {
    int c_outlier_idx = dis2(gen);
    assert(c_outlier_idx < expected_outlier_mask.size());
    if(expected_outlier_mask[c_outlier_idx])
    {
      i--;
      continue;
    }
    expected_outlier_mask[c_outlier_idx] = true;
    Eigen::Matrix<double, 3, 1> Tt;
    for(int tt = 0; tt < 3; tt++)
    {
      if(disC(gen) <= 0.5)
        Tt(tt) = dis3_lb(gen);
      else
        Tt(tt) = dis3_ub(gen);
    }
    tgt.col(c_outlier_idx) += Tt; // random translation
  }
}

int main(int argv, char* argc[]) {
  double sumSE = 0, sumRE = 0, sumTE = 0, sumTime = 0, sumRMSE;
  int totalData = argv - 1;
  for(int t = 1; t < argv; t++)
  {
    // Load the .ply file
    teaser::PLYReader reader;
    teaser::PointCloud src_cloud;
    auto status = reader.read(std::string(argc[t]), src_cloud);
    int N = src_cloud.size();

    // Convert the point cloud to Eigen
    Eigen::Matrix<double, 3, Eigen::Dynamic> src(3, N);
    for (size_t i = 0; i < N; ++i) {
      src.col(i) << src_cloud[i].x, src_cloud[i].y, src_cloud[i].z;
    }

    // Homogeneous coordinates
    Eigen::Matrix<double, 4, Eigen::Dynamic> src_h;
    src_h.resize(4, src.cols());
    src_h.topRows(3) = src;
    src_h.bottomRows(1) = Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(N);

    srand( time(NULL) );
    //test mutiple + excel
    double scale_error_sum = 0, scale_error_sum_of_square = 0;
    double angle_error_sum = 0, angle_error_sum_of_square = 0;
    double trans_error_sum = 0, trans_error_sum_of_square = 0;
    double time_sum = 0, time_sum_of_square = 0;
    double RMSE_sum = 0, RMSE_sum_of_square = 0;
    int testTime = 100;
    for(int i = 0; i < testTime; i++)
    {      
      //random
      double testScale = 1;
      if (unknownScale)
        testScale = 1 + 4 * ((1.0 - 0.0) * rand() / (RAND_MAX + 1.0) + 0.0);
      Eigen::Matrix4d T;
      T << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1;
      Eigen::Matrix<double, 3, Eigen::Dynamic> axis = Eigen::Matrix<double, 3, Eigen::Dynamic>::Random(3, 1);
      axis = axis/axis.norm();
      double randomAngle = (3.1416 - 0.0) * rand() / (RAND_MAX + 1.0) + 0.0;
      Eigen::Matrix<double, 3, Eigen::Dynamic> aa = randomAngle * axis;
      Eigen::Matrix<double, 3, 3> R;
      R << 1, 0, 0,
          0, 1, 0,
          0, 0, 1;
      if(aa.norm() >= 2e-16)
      {
        Eigen::Matrix<double, 3, Eigen::Dynamic> K = aa / aa.norm();
        Eigen::Matrix<double, 3, 3> K1;
        K1 << 0, -K(2), K(1),
              K(2), 0, -K(0),
              -K(1), K(0), 0;
        R = R + sin(aa.norm()) * K1 + (1 - cos(aa.norm())) * K1 * K1;
      }
      T.topLeftCorner(3, 3) = R;
      Eigen::Matrix<double, 3, 1> Tt;
      Tt(0) = (1.0 - 0.0) * rand() / (RAND_MAX + 1.0) + 0.0 - 0.5;
      Tt(1) = (1.0 - 0.0) * rand() / (RAND_MAX + 1.0) + 0.0 - 0.5;
      Tt(2) = (1.0 - 0.0) * rand() / (RAND_MAX + 1.0) + 0.0 - 0.5;
      Tt = 3 * ((1.0 - 0.0) * rand() / (RAND_MAX + 1.0) + 0.0) * Tt / Tt.norm();
      T.topRightCorner(3, 1) = Tt;

      // Apply transformation
      Eigen::Matrix<double, 4, Eigen::Dynamic> tgt_h = testScale * T * src_h;
      Eigen::Matrix<double, 3, Eigen::Dynamic> tgt = tgt_h.topRows(3);
      Eigen::Matrix<double, 3, Eigen::Dynamic> src_gt = tgt_h.topRows(3);

      // Add some noise & outliers
      addNoiseAndOutliers(tgt);

      // Run our modified TEASER++ registration
      // Prepare solver parameters
      teaser::RobustRegistrationSolver::Params params;
      params.noise_bound = NOISE_BOUND;
      params.cbar2 = 1;
      params.estimate_scaling = unknownScale;
      params.rotation_max_iterations = 100;
      params.rotation_gnc_factor = 1.4;
      params.rotation_estimation_algorithm =
          teaser::RobustRegistrationSolver::ROTATION_ESTIMATION_ALGORITHM::GNC_TLS;
      params.rotation_cost_threshold = 0.005;

      // Solve with our modified TEASER++ method with the propose C-RANSAC structure
      teaser::RobustRegistrationSolver solver(params);
      std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
      solver.solve(src, tgt);
      std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

      auto solution = solver.getSolution();
      
      double angleError = getAngularError(T.topLeftCorner(3, 3), solution.rotation);
      scale_error_sum += testScale - solution.scale;
      scale_error_sum_of_square += pow(testScale - solution.scale, 2);
      angle_error_sum += angleError;
      angle_error_sum_of_square += pow(angleError, 2);
      trans_error_sum += (T.topRightCorner(3, 1) - solution.translation).norm();
      trans_error_sum_of_square += pow((T.topRightCorner(3, 1) - solution.translation).norm(), 2);
      time_sum += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
      time_sum_of_square += pow(std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0, 2);
      // Homogeneous coordinates
      Eigen::Matrix4d TRANSFORM;
      TRANSFORM <<  1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1;
      TRANSFORM.topLeftCorner(3, 3) = solution.rotation;
      TRANSFORM.topRightCorner(3, 1) = solution.translation;
      Eigen::Matrix<double, 3, Eigen::Dynamic> src_solve = (solution.scale * TRANSFORM * src_h).topRows(3);
      double residualSum = 0.0;
      for(int m = 0; m < src_solve.cols(); m++)
      {
        residualSum += pow((src_gt.col(m) - src_solve.col(m)).norm(), 2);
      }
      residualSum = sqrt(residualSum / src_solve.cols());
      RMSE_sum += residualSum;
      RMSE_sum_of_square += pow(residualSum, 2);

      // Compare results
      std::cout << "\n\n";
      std::cout << "=====================================" << std::endl;
      std::cout << "          C-RANSAC Results           " << std::endl;
      std::cout << "=====================================" << std::endl;
      std::cout << "Expected scale: " << std::endl;
      std::cout << testScale << std::endl;
      std::cout << "Estimated scale: " << std::endl;
      std::cout << solution.scale << std::endl;
      std::cout << "Error (): " << testScale - solution.scale
                << std::endl;
      std::cout << std::endl;
      std::cout << "Expected rotation: " << std::endl;
      std::cout << T.topLeftCorner(3, 3) << std::endl;
      std::cout << "Estimated rotation: " << std::endl;
      std::cout << solution.rotation << std::endl;
      std::cout << "Error (deg): " << angleError
                << std::endl;
      std::cout << std::endl;
      std::cout << "Expected translation: " << std::endl;
      std::cout << T.topRightCorner(3, 1) << std::endl;
      std::cout << "Estimated translation: " << std::endl;
      std::cout << solution.translation << std::endl;
      std::cout << "Error (m): " << (T.topRightCorner(3, 1) - solution.translation).norm() << std::endl;
      std::cout << std::endl;
      std::cout << "RMSE: " << residualSum << std::endl;
      std::cout << std::endl;
      std::cout << "Number of correspondences: " << N << std::endl;
      std::cout << "Number of outliers: " << tgt.cols() * N_OUTLIERS_RATE << std::endl;
      std::cout << "Time taken (s): "
                << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() /
                      1000000.0
                << std::endl;
      std::cout << std::string(argc[t]) << " " << i << std::endl;
      sleep(3);
    }
    double scale_error_standard_deviation = sqrt(scale_error_sum_of_square / testTime - pow(scale_error_sum / testTime, 2)); //calculate standard deviation
    double angle_error_standard_deviation = sqrt(angle_error_sum_of_square / testTime - pow(angle_error_sum / testTime, 2)); //calculate standard deviation
    double trans_error_standard_deviation = sqrt(trans_error_sum_of_square / testTime - pow(trans_error_sum / testTime, 2)); //calculate standard deviation
    double time_standard_deviation = sqrt(time_sum_of_square / testTime - pow(time_sum / testTime, 2)); //calculate standard deviation
    double RMSE_standard_deviation = sqrt(RMSE_sum_of_square / testTime - pow(RMSE_sum / testTime, 2)); //calculate standard deviation
    
    std::cout << "=====================================" << std::endl;
    std::cout << "                Final                " << std::endl;
    std::cout << "=====================================" << std::endl;
    std::cout << "Scale_Error:, (esti. / grou.)" << std::endl;
    std::cout << scale_error_sum / testTime + 2 * scale_error_standard_deviation << std::endl;
    std::cout << scale_error_sum / testTime + 1 * scale_error_standard_deviation << std::endl;
    std::cout << scale_error_sum / testTime << std::endl;
    std::cout << scale_error_sum / testTime - 1 * scale_error_standard_deviation << std::endl;
    std::cout << scale_error_sum / testTime - 2 * scale_error_standard_deviation << std::endl;
    std::cout << std::endl << "Angle_Error:, (degree)" << std::endl;
    std::cout << angle_error_sum / testTime + 2 * angle_error_standard_deviation << std::endl;
    std::cout << angle_error_sum / testTime + 1 * angle_error_standard_deviation << std::endl;
    std::cout << angle_error_sum / testTime << std::endl;
    std::cout << angle_error_sum / testTime - 1 * angle_error_standard_deviation << std::endl;
    std::cout << angle_error_sum / testTime - 2 * angle_error_standard_deviation << std::endl;
    std::cout << std::endl << "Trans_Error:, (meter)" << std::endl;
    std::cout << trans_error_sum / testTime + 2 * trans_error_standard_deviation << std::endl;
    std::cout << trans_error_sum / testTime + 1 * trans_error_standard_deviation << std::endl;
    std::cout << trans_error_sum / testTime << std::endl;
    std::cout << trans_error_sum / testTime - 1 * trans_error_standard_deviation << std::endl;
    std::cout << trans_error_sum / testTime - 2 * trans_error_standard_deviation << std::endl;
    std::cout << std::endl << "Time:, (second)" << std::endl;
    std::cout << time_sum / testTime + 2 * time_standard_deviation << std::endl;
    std::cout << time_sum / testTime + 1 * time_standard_deviation << std::endl;
    std::cout << time_sum / testTime << std::endl;
    std::cout << time_sum / testTime - 1 * time_standard_deviation << std::endl;
    std::cout << time_sum / testTime - 2 * time_standard_deviation << std::endl;
    std::cout << std::endl << "RMSE" << std::endl;
    std::cout << RMSE_sum / testTime + 2 * RMSE_standard_deviation << std::endl;
    std::cout << RMSE_sum / testTime + 1 * RMSE_standard_deviation << std::endl;
    std::cout << RMSE_sum / testTime << std::endl;
    std::cout << RMSE_sum / testTime - 1 * RMSE_standard_deviation << std::endl;
    std::cout << RMSE_sum / testTime - 2 * RMSE_standard_deviation << std::endl;
    //excel
    std::string s(argc[t]); s.resize(s.length() - 4);
    std::ofstream oFile;
    oFile.open(s + ".csv", std::ios::out | std::ios::trunc);
    oFile << "Scale_Error:, (esti. / grou.)" << std::endl; 
    oFile << "m + 2v," << scale_error_sum / testTime + 2 * scale_error_standard_deviation << std::endl;
    oFile << "m + 1v," << scale_error_sum / testTime + 1 * scale_error_standard_deviation << std::endl;
    oFile << "m," << scale_error_sum / testTime << std::endl;
    oFile << "m - 1v," << scale_error_sum / testTime - 1 * scale_error_standard_deviation << std::endl;
    oFile << "m - 2v," << scale_error_sum / testTime - 2 * scale_error_standard_deviation << std::endl;
    oFile << std::endl << "Angle_Error:, (degree)" << std::endl;
    oFile << "m + 2v,"<< angle_error_sum / testTime + 2 * angle_error_standard_deviation << std::endl;
    oFile << "m + 1v,"<< angle_error_sum / testTime + 1 * angle_error_standard_deviation << std::endl;
    oFile << "m,"<< angle_error_sum / testTime << std::endl;
    oFile << "m - 1v,"<< angle_error_sum / testTime - 1 * angle_error_standard_deviation << std::endl;
    oFile << "m - 2v,"<< angle_error_sum / testTime - 2 * angle_error_standard_deviation << std::endl;
    oFile << std::endl << "Trans_Error:, (meter)" << std::endl;
    oFile << "m + 2v,"<< trans_error_sum / testTime + 2 * trans_error_standard_deviation << std::endl;
    oFile << "m + 1v,"<< trans_error_sum / testTime + 1 * trans_error_standard_deviation << std::endl;
    oFile << "m,"<< trans_error_sum / testTime << std::endl;
    oFile << "m - 1v,"<< trans_error_sum / testTime - 1 * trans_error_standard_deviation << std::endl;
    oFile << "m - 2v,"<< trans_error_sum / testTime - 2 * trans_error_standard_deviation << std::endl;
    oFile << std::endl << "Time:, (second)" << std::endl;
    oFile << "m + 2v,"<< time_sum / testTime + 2 * time_standard_deviation << std::endl;
    oFile << "m + 1v,"<< time_sum / testTime + 1 * time_standard_deviation << std::endl;
    oFile << "m,"<< time_sum / testTime << std::endl;
    oFile << "m - 1v,"<< time_sum / testTime - 1 * time_standard_deviation << std::endl;
    oFile << "m - 2v,"<< time_sum / testTime - 2 * time_standard_deviation << std::endl;
    oFile << std::endl << "RMSE" << std::endl;
    oFile << "m + 2v,"<< RMSE_sum / testTime + 2 * RMSE_standard_deviation << std::endl;
    oFile << "m + 1v,"<< RMSE_sum / testTime + 1 * RMSE_standard_deviation << std::endl;
    oFile << "m,"<< RMSE_sum / testTime << std::endl;
    oFile << "m - 1v,"<< RMSE_sum / testTime - 1 * RMSE_standard_deviation << std::endl;
    oFile << "m - 2v,"<< RMSE_sum / testTime - 2 * RMSE_standard_deviation << std::endl;
    oFile.close();
    std::cout << std::string(argc[t]) << " End." << std::endl;

    sumSE += scale_error_sum / testTime;
    sumRE += angle_error_sum / testTime;
    sumTE += trans_error_sum / testTime;
    sumTime += time_sum / testTime;
    sumRMSE += RMSE_sum / testTime;     
    sleep(3);
  }
  
  std::ofstream oFile;
  oFile.open("Average.csv", std::ios::out | std::ios::trunc);
  oFile << "Scale_Error:, (esti. / grou.)" << std::endl;
  oFile << "," << sumSE / totalData << std::endl;
  oFile << "Angle_Error:, (degree)" << std::endl;
  oFile << ","<< sumRE / totalData << std::endl;
  oFile << "Trans_Error:, (meter)" << std::endl;
  oFile << ","<< sumTE / totalData << std::endl;
  oFile << "Time:, (second)" << std::endl;
  oFile << ","<< sumTime / totalData << std::endl;
  oFile << "RMSE:" << std::endl;
  oFile << ","<< sumRMSE / totalData << std::endl;
  oFile.close();  
}
