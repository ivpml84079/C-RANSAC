// An example showing TEASER++ registration with FPFH features with the Stanford bunny model
#include <chrono>
#include <iostream>
#include <random>

#include <Eigen/Core>

#include <teaser/ply_io.h>
#include <teaser/registration.h>
#include <teaser/matcher.h>

#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/octree/octree.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/features/normal_3d.h>
#include <pcl/search/kdtree.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/filter.h>
#include <pcl/registration/icp.h>
#include <pcl/registration/ndt.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/common/transformation_from_correspondences.h>
#include <pcl/kdtree/kdtree_flann.h>


// Macro constants for generating noise and outliers
#define NOISE_BOUND 0.05
#define N_OUTLIERS 1700
#define OUTLIER_TRANSLATION_LB 5
#define OUTLIER_TRANSLATION_UB 10

inline double getAngularError(Eigen::Matrix3d R_exp, Eigen::Matrix3d R_est) {
  return std::abs(std::acos(fmin(fmax(((R_exp.transpose() * R_est).trace() - 1) / 2, -1.0), 1.0)));
}

void addNoiseAndOutliers(Eigen::Matrix<double, 3, Eigen::Dynamic>& tgt) {
  // Add uniform noise
  Eigen::Matrix<double, 3, Eigen::Dynamic> noise =
      Eigen::Matrix<double, 3, Eigen::Dynamic>::Random(3, tgt.cols()) * NOISE_BOUND;
  NOISE_BOUND / 2;
  tgt = tgt + noise;

  // Add outliers
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis2(0, tgt.cols() - 1); // pos of outliers
  std::uniform_int_distribution<> dis3(OUTLIER_TRANSLATION_LB,
                                       OUTLIER_TRANSLATION_UB); // random translation
  std::vector<bool> expected_outlier_mask(tgt.cols(), false);
  for (int i = 0; i < N_OUTLIERS; ++i) {
    int c_outlier_idx = dis2(gen);
    assert(c_outlier_idx < expected_outlier_mask.size());
    expected_outlier_mask[c_outlier_idx] = true;
    tgt.col(c_outlier_idx).array() += dis3(gen); // random translation
  }
}
void compute_fpfh(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, double normal_search_radius, double fpfh_search_radius, pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfh_out){
  pcl::PointCloud<pcl::Normal>::Ptr normal(new pcl::PointCloud<pcl::Normal>());
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
	ne.setInputCloud(cloud);
	ne.setSearchMethod(tree);
	ne.setRadiusSearch(normal_search_radius);
	ne.compute(*normal);

	//compute fpfh using normals
	pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_est;
	fpfh_est.setInputCloud(cloud);
	fpfh_est.setInputNormals(normal);
	fpfh_est.setSearchMethod(tree);
	fpfh_est.setRadiusSearch(fpfh_search_radius);
	fpfh_est.compute(*fpfh_out);
}
void correspondenceSearching(pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs, pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfht, pcl::Correspondences & corr, int max_corr, std::vector<int>& corr_NOs, std::vector<int>& corr_NOt)
{
	int n = std::min(max_corr, (int)fpfht->size()); //maximum number of correspondences to find for each source point
	corr.clear();
	corr_NOs.assign(fpfhs->size(), 0);
	corr_NOt.assign(fpfht->size(), 0);
	// Use a KdTree to search for the nearest matches in feature space
	pcl::KdTreeFLANN<pcl::FPFHSignature33> treeS;
	treeS.setInputCloud(fpfhs);
	pcl::KdTreeFLANN<pcl::FPFHSignature33> treeT;
	treeT.setInputCloud(fpfht);
	for (size_t i = 0; i < fpfhs->size(); i++) {
		std::vector<int> corrIdxTmp(n);
		std::vector<float> corrDisTmp(n);
		//find the best n matches in target fpfh
		treeT.nearestKSearch(*fpfhs, i, n, corrIdxTmp, corrDisTmp);
		for (size_t j = 0; j < corrIdxTmp.size(); j++) {
			bool removeFlag = true;
			int searchIdx = corrIdxTmp[j];
			std::vector<int> corrIdxTmpT(n);
			std::vector<float> corrDisTmpT(n);
			treeS.nearestKSearch(*fpfht, searchIdx, n, corrIdxTmpT, corrDisTmpT);
			for (size_t k = 0; k < n; k++) {
				if (corrIdxTmpT.data()[k] == i) {
					removeFlag = false;
					break;
				}
			}
			if (removeFlag == false) {
				pcl::Correspondence corrTabTmp;
				corrTabTmp.index_query = i;
				corrTabTmp.index_match = corrIdxTmp[j];
				corrTabTmp.distance = corrDisTmp[j];
				corr.push_back(corrTabTmp);
				corr_NOs[i]++;
				corr_NOt[corrIdxTmp[j]]++;
			}
		}
	}
}
void search_corres_fpfh(pcl::PointCloud<pcl::PointXYZ>::Ptr src, pcl::PointCloud<pcl::PointXYZ>::Ptr tar, std::vector<std::pair<int, int>> &corr, double normal_search_radius, double fpfh_search_radius){
  pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhS(new pcl::PointCloud<pcl::FPFHSignature33>());
	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhT(new pcl::PointCloud<pcl::FPFHSignature33>());
  compute_fpfh(src, normal_search_radius, fpfh_search_radius, fpfhS);
  compute_fpfh(tar, normal_search_radius, fpfh_search_radius, fpfhT);

  int max_corr = 5; // 只找最高分的 5 個
  std::vector<int> corr_NOS, corr_NOT;
  pcl::CorrespondencesPtr corr_temp(new pcl::Correspondences);
	correspondenceSearching(fpfhS, fpfhT, *corr_temp, max_corr, corr_NOS, corr_NOT);

  int N = corr_NOS.size();
  corr.resize(N);
  for(int i=0; i<N; i++){
    corr[i] = {(*corr_temp)[i].index_query, (*corr_temp)[i].index_match};
  }

  return;
}

int main() {
  // Load the .ply file
  // teaser::PLYReader reader;
  // teaser::PointCloud src_cloud;
  // auto status = reader.read("/mnt/e/Chris/DataSet/SS_exp/RESSO/figure_7a/part0.ply", src_cloud);
  // auto status = reader.read("./example_data/bun_zipper_res3.ply", src_cloud);
  // int N = src_cloud.size();

  // Convert the point cloud to Eigen
  // Eigen::Matrix<double, 3, Eigen::Dynamic> src(3, N);
  // for (size_t i = 0; i < N; ++i) {
  //   src.col(i) << src_cloud[i].x, src_cloud[i].y, src_cloud[i].z;
  // }

  // std::string fnameS = "/mnt/e/Chris/DataSet/SS_exp/office(re)/s1.ply";
  // std::string fnameT = "/mnt/e/Chris/DataSet/SS_exp/office(re)/s5.ply";
  // std::string fnameS = "/mnt/e/Chris/DataSet/SS_exp/RESSO/figure_7a/part1.ply";
  // std::string fnameT = "/mnt/e/Chris/DataSet/SS_exp/RESSO/figure_7a/part0.ply";
  std::string fnameS = "/mnt/d/2022/[PAPER]/3D_point_cloud/[DATASET]/ETH_auto_reg/trees/s1.ply";
  std::string fnameT = "/mnt/d/2022/[PAPER]/3D_point_cloud/[DATASET]/ETH_auto_reg/trees/s2.ply";
  // std::string fnameS = "/mnt/e/Chris/DataSet/SS_exp/RESSO/figure_6i/part4.ply";
  // std::string fnameT = "/mnt/e/Chris/DataSet/SS_exp/RESSO/figure_6i/part0.ply";

  // std::string fnameS = "./example_data/bun_zipper_res3.ply";
  // std::string fnameT = "./example_data/bun_zipper_res3.ply";

  pcl::PointCloud<pcl::PointXYZ>::Ptr source(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr target(new pcl::PointCloud<pcl::PointXYZ>);

  if (pcl::io::loadPLYFile<pcl::PointXYZ>(fnameS, *source) == -1)
	{
		PCL_ERROR("Couldn't read file \n");
		return 0;
	}

	if (pcl::io::loadPLYFile<pcl::PointXYZ>(fnameT, *target) == -1)
	{
		PCL_ERROR("Couldn't read file \n");
		return 0;
	}
  std::cout << "reading complete\n";

  // Preprocessing (downsampling)
  float LeafSize = 0.2;
  pcl::VoxelGrid<pcl::PointXYZ> voxel_grid_source;
  voxel_grid_source.setLeafSize(LeafSize, LeafSize, LeafSize);
  voxel_grid_source.setInputCloud(source);
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src(new pcl::PointCloud<pcl::PointXYZ>);
  voxel_grid_source.filter(*cloud_src);
  
  pcl::VoxelGrid<pcl::PointXYZ> voxel_grid_target;
  voxel_grid_source.setLeafSize(LeafSize, LeafSize, LeafSize);
  voxel_grid_source.setInputCloud(target);
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tar(new pcl::PointCloud<pcl::PointXYZ>);
  voxel_grid_source.filter(*cloud_tar);
  std::cout << "downsampling complete\n";

  // Conver the pcl point cloud to teaser point cloud
  int N_src = cloud_src->points.size();
  int N_tar = cloud_tar->points.size();
  teaser::PointCloud src_cloud;
  teaser::PointCloud tgt_cloud;
  for(int i=0; i<N_src; i++){
    pcl::PointXYZ temp_point = cloud_src->points[i];
    src_cloud.push_back({temp_point.x, temp_point.y, temp_point.z});
  }
  for(int i=0; i<N_tar; i++){
    pcl::PointXYZ temp_point = cloud_tar->points[i];
    tgt_cloud.push_back({temp_point.x, temp_point.y, temp_point.z});
  }
  std::cout << "convert complete\n";
  std::cout << "src cloud : " << N_src << std::endl;
  std::cout << "tar cloud : " << N_tar << std::endl;


  // Convert the point cloud to Eigen

  // Eigen::Matrix<double, 3, Eigen::Dynamic> src(3, N_src);
  // for (size_t i = 0; i < N_src; ++i) {
  //   src.col(i) << src_cloud->points[i].x, src_cloud->points[i].y, src_cloud->points[i].z;
  // }
  // Eigen::Matrix<double, 3, Eigen::Dynamic> tgt(3, N_tar);
  // for (size_t i = 0; i < N_tar; ++i) {
  //   tgt.col(i) << tgt_cloud->points[i].x, tgt_cloud->points[i].y, tgt_cloud->points[i].z;
  // }



  // std::cout << "read point cloud with " << N_src << " points\n";
  // std::cout << "read point cloud with " << N_tar << " points\n";

  // Homogeneous coordinates
  // Eigen::Matrix<double, 4, Eigen::Dynamic> src_h;
  // src_h.resize(4, src.cols());
  // src_h.topRows(3) = src;
  // src_h.bottomRows(1) = Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(N);

  // Apply an arbitrary SE(3) transformation
  Eigen::Matrix4d T;
  // clang-format off
  
  // GT 放在T裡面
  // GT for Figure 7a.
  // T << 1, 0, 0, 0,
  //       0, 1, 0 ,0,
  //       0, 0, 1, 0,
  //       0, 0, 0, 1;
  // GT for Office s1-s5
  // T << -0.5769168846,-0.8156505245,0.0433719971,-3.6717774167,
  //       0.8167742979,-0.5765268750,0.0222824749,3.9496681535,
  //       0.0068304096,0.0482802685,0.9988104731,0.0620916024,
  //       0.0000000000,0.0000000000,0.0000000000,1.0000000000;
  // GT for Heritage Building
  T << -0.71860257, -0.69542063, 0.00069803, -11.10224941,
        0.69541843, -0.7186023, -0.00198637, -0.8739701127,
        0.00188297, -0.00094199, 0.99999778, 0.5416248644,
        0, 0, 0, 1;
  // T <<  0.227703, 0.973731, -0.000262, -22.620619, 
	// 		 -0.973717, 0.227701, 0.005246, -3.725392, 
  //       0.005167, -0.000940, 0.999986, 0.225658,
  //       0, 0, 0, 1;
  // GT for Resso6i 0-4
	// T <<  0.979955,    0.195595,  -0.0378434, -2.09979,
  // 				-0.195148,    0.980656,   0.0152146, 1.76114,
  // 				0.0400873, -0.00752457,    0.999169, -0.0496839,
  //         0, 0, 0, 1;

  // clang-format on

  // Apply transformation
  // Eigen::Matrix<double, 4, Eigen::Dynamic> tgt_h = T * src_h;
  // Eigen::Matrix<double, 3, Eigen::Dynamic> tgt = tgt_h.topRows(3);

  // Add some noise & outliers
  // addNoiseAndOutliers(tgt);

  // Convert to teaser point cloud
  // teaser::PointCloud tgt_cloud;
  // for (size_t i = 0; i < tgt.cols(); ++i) {
  //   tgt_cloud.push_back({static_cast<float>(tgt(0, i)), static_cast<float>(tgt(1, i)),
  //                        static_cast<float>(tgt(2, i))});
  // }

  // Compute FPFH
  teaser::FPFHEstimation fpfh;
  auto obj_descriptors = fpfh.computeFPFHFeatures(src_cloud, 3.0, 3.0);
  auto scene_descriptors = fpfh.computeFPFHFeatures(tgt_cloud, 3.0, 3.0);

  teaser::Matcher matcher;
  clock_t f_s, f_t;
  f_s = clock();
  auto correspondences = matcher.calculateCorrespondences(
      src_cloud, tgt_cloud, *obj_descriptors, *scene_descriptors, false, true, false, 0.95);
  f_t = clock();
  std::cout << "FPFH cost " << double(f_t-f_s) / CLOCKS_PER_SEC << std::endl;


  // Compute FPFH(self)
  // std::vector<std::pair<int, int>> correspondences;
  // search_corres_fpfh(src_cloud, tgt_cloud, correspondences, 0.02, 0.04);

  // Run TEASER++ registration
  // Prepare solver parameters
  teaser::RobustRegistrationSolver::Params params;
  params.noise_bound = NOISE_BOUND;
  params.cbar2 = 1;
  params.estimate_scaling = false;
  params.rotation_max_iterations = 100;
  params.rotation_gnc_factor = 1.4;
  params.rotation_estimation_algorithm =
      teaser::RobustRegistrationSolver::ROTATION_ESTIMATION_ALGORITHM::GNC_TLS;
  params.rotation_cost_threshold = 0.005;

  float mean_r = 0.0;
  float mean_t = 0.0;
  float mean_time = 0.0;
  for(int test_time=0; test_time<100; test_time++){
    // Solve with TEASER++
    teaser::RobustRegistrationSolver solver(params);
    // std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    clock_t t1, t2;
    t1 = clock();
    solver.solve(src_cloud, tgt_cloud, correspondences);
    t2 = clock();
    // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    auto solution = solver.getSolution();

    // std::cout << "solution : " << solution << std::endl;

    // Compare results
    // std::cout << "=====================================" << std::endl;
    // std::cout << "          TEASER++ Results           " << std::endl;
    // std::cout << "=====================================" << std::endl;
    // // std::cout << "Expected rotation: " << std::endl;
    // std::cout << T.topLeftCorner(3, 3) << std::endl;
    // std::cout << "Estimated rotation: " << std::endl;
    // std::cout << solution.rotation << std::endl;
    // std::cout << "Error (deg): " << getAngularError(T.topLeftCorner(3, 3), solution.rotation)
    //           << std::endl;
    // std::cout << std::endl;
    // std::cout << "Expected translation: " << std::endl;
    // std::cout << T.topRightCorner(3, 1) << std::endl;
    // std::cout << "Estimated translation: " << std::endl;
    // std::cout << solution.translation << std::endl;
    // std::cout << "Error (m): " << (T.topRightCorner(3, 1) - solution.translation).norm() << std::endl;
    // std::cout << std::endl;
    // std::cout << "Number of correspondences: " << N << std::endl;
    // std::cout << "Number of outliers: " << N_OUTLIERS << std::endl;
    // std::cout << "Time taken (s): "
    //           << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() /
    //                 1000000.0
    //           << std::endl;
    mean_r += getAngularError(T.topLeftCorner(3, 3), solution.rotation);
    mean_t += (T.topRightCorner(3, 1) - solution.translation).norm();
    // mean_time += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
    mean_time += double(t2-t1) / CLOCKS_PER_SEC;
  }
  std::cout << "r : " << mean_r / 100 << std::endl;
  std::cout << "t : " << mean_t / 100 << std::endl;
  std::cout << "time : " << mean_time / 100 << std::endl;
}