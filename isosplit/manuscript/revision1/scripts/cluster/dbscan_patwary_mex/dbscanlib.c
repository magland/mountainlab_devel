// C callable interface to DBSCAN by Patwary
// Barnett 4/17/15

#include "dbscan.h"
#include "utils.h"
#include "kdtree2.hpp"

void dbscanwrapper(float** X, int dims, int N, int minPts, double eps,
	       int* K, int* L)
/*

  L - labels output (zero if not clus)
  K - # found
*/
{

  NWUClustering::ClusteringAlgo dbs;
  
  dbs.set_dbscan_params(eps, minPts);

  // hack dbs.read_file
  
  dbs.build_kdtree();

  run_dbscan_algo_uf(dbs);

  // hack dbs.writeClusters_uf
  
  return 0;
}
