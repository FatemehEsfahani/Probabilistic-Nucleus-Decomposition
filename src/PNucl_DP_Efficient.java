//This version uses dynamic programming. Also, it considers all the triangles with 
//existence probability greater than eta.

import it.unimi.dsi.webgraph.ImmutableGraph;
import it.unimi.dsi.webgraph.labelling.ArcLabelledImmutableGraph;
import it.unimi.dsi.webgraph.labelling.Label;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;

import java.util.List;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.stream.IntStream;

public class PNucl_DP_Efficient {

	String basename; // uncertain graph
	int precision;
	double eta;
	int num_tria_alive;
	
	double[] triangle_prob;

	ArcLabelledImmutableGraph G;
	int n; // vertices of G
	int m; // edges of G

	BigInteger t_total_cnt = BigInteger.ZERO;
	int maxdeg;

	double epsilon0 = 1.0E-16;
	double epsilon1 = 1.0 - epsilon0;

	int[] _4kDegree; // stores 4kdegree for each triangle

	BitSet processed; // records if a triangle has been processed or not
	int[] sortedTrinagle; // stores triangles in ascending order of their 4kdegrees
	int[] TrianglePos; // stores position of each triangle in sortedTrinagle array

	// for each triangle (u,v,w), where u<v<w:
	int[] first_vert; // first_vert stores u
	int[] second_vert; // second_vert stores v
	int[] third_vert; // third_vert stores
	int[] offset_uv; // offset array for triangles

	final int[] edgeTail;
	final int[] offset; // offset array for edges
	double[] edge_prob;

	public PNucl_DP_Efficient(String basename, int precision, double eta) throws Exception {

		this.basename = basename;
		this.precision = precision;
		this.eta = eta;

		System.out.println("loading the graph...");
		long startTime1 = System.currentTimeMillis();
		G = ArcLabelledImmutableGraph.load(basename);
		n = G.numNodes();
		m = (int) (G.numArcs() / 2);
		long endTime1 = System.currentTimeMillis();
		System.out.println("m " + m + " n " + n);
		System.out.println("Time elapsed (sec) for loading the graph = " + (endTime1 - startTime1) / 1000.0);

		offset_uv = new int[m + 1]; // starting index of triangles with edge (u,v) where u<v

		edgeTail = new int[m]; // store tail for each edge
		offset = new int[n + 1]; // edge offset for edgeHead array
		edge_prob = new double[m];

		System.out.println("Start storing edges in edgeTail.");
		long startTime2 = System.currentTimeMillis();
		Populate_edges();
		long endTime2 = System.currentTimeMillis();
		System.out.println("Time elapsed (sec) for storing edges = " + (endTime2 - startTime2) / 1000.0);

		
		EnumerateTriangles(); // In this method the value of "num_tria_alive" is set.

		first_vert = new int[num_tria_alive];
		second_vert = new int[num_tria_alive];
		third_vert = new int[num_tria_alive];

		_4kDegree = new int[num_tria_alive];
		triangle_prob = new double[num_tria_alive];

		processed = new BitSet(num_tria_alive); // no triangle has been processed yet

		sortedTrinagle = new int[num_tria_alive];
		TrianglePos = new int[num_tria_alive];

		long startTime3 = System.currentTimeMillis();
		Populate_triangles();
		long endTime3 = System.currentTimeMillis();
		System.out.println("Time elapsed (sec) for storing trinagles = " + (endTime3 - startTime3) / 1000.0);

		// This pass is very fast as it does not read the neighbors of vertices.
		maxdeg = 0;
		for (int v = 0; v < n; v++) {
			int v_deg = G.outdegree(v);
			if (v_deg > maxdeg)
				maxdeg = v_deg;
		}

		System.out.println("n=" + n + ", m=" + m + ", maxdeg=" + maxdeg);
		
		
	}

	// outputs the number of triangles greater than eta (threshold)
	public void EnumerateTriangles() throws Exception {

		System.out.println("start enumearting trinagles...");
		t_total_cnt = IntStream.range(0, n).map(u -> {
			if (u % 1000 == 0)
				System.out.println(u);

			ImmutableGraph H = G.copy();
			int[] u_neighbors = H.successorArray(u);
			int u_deg = H.outdegree(u);

			int t_cnt = IntStream.range(0, u_deg).map(i -> {
				int v = u_neighbors[i];
				if (u < v) { // to avoid double counting - for non-optimized graph

					int[] v_neighbors = H.successorArray(v);
					int v_deg = H.outdegree(v);
					List<Integer> intersection = intersection(u, v, u_neighbors, v_neighbors, u_deg, v_deg);
					int index = intersection.size();

					return index; // this is just the count of triangles
				} else {
					return 0;
				}
			}).reduce(0, (a, b) -> a + b);
			return t_cnt;
		}).mapToObj(BigInteger::valueOf).reduce(BigInteger.ZERO, BigInteger::add);

		System.out.println("Number of Triangles with high existence probability: " + t_total_cnt);
		this.num_tria_alive = t_total_cnt.intValue();
	}

	// stores edges of the graph
	public void Populate_edges() {
		int startPos = 0;
		for (int u = 0; u < n; u++) {
			offset[u] = startPos;

			int[] u_successor = G.successorArray(u);
			Label[] u_label = G.labelArray(u);

			// Arrays.sort(u_successor);

			int count = 0;
			for (int j = 0; j < u_successor.length; j++) {
				int v = u_successor[j];
				if (u < v) {
					edgeTail[startPos + count] = v;
					edge_prob[startPos + count] = u_label[j].getLong() * (1.0 / Math.pow(10, precision));
					//System.out.println(edge_prob[startPos + count]);
					count++;
				}
			}
			startPos += count;
		}
		offset[n] = offset[n - 1];
		
		 //double sum = 0.0;
		//double sum2 = 0.0;
		
		//int count =0;
		//for (int i=0;i<m;i++) {
		//	sum2 += edge_prob[i];
			//if (edge_prob[i] > 0.5 ) {
			//	sum += edge_prob[i];
				//count++;
			//}
			
		//}
		
		//double ave = sum/count;
		//System.out.println( ave + " " + (100*count/(double) m));
		//System.out.println("average edge probability is: " + sum2/m);
	}

	// stores triangles of the graph
	public void Populate_triangles() {
		int startPos = 0;
		int count_triangles_ddetected_so_far = 0;

		for (int e = 0; e < m; e++) {
			int u = findHead(e);
			int v = edgeTail[e];
			offset_uv[e] = startPos;

			int[] u_neighbors = G.successorArray(u);
			int[] v_neighbors = G.successorArray(v);

			int u_deg = G.outdegree(u);
			int v_deg = G.outdegree(v);

			int index = intersection2(u, v, u_neighbors, v_neighbors, u_deg, v_deg, count_triangles_ddetected_so_far);
			// for each edge (u,v) index is the number of triangles containing edge (u,v)

			count_triangles_ddetected_so_far += index;
			startPos += index;
		}

		offset_uv[m] = offset_uv[m - 1];

	}

	// computes initial degree for each triangle
	public void computeTriangles_degree() throws Exception {

		System.out.println("start computing dgeree for each triangle");

		long startTime = System.currentTimeMillis();

		t_total_cnt = IntStream.range(0, n).parallel().map(u -> {

			ArcLabelledImmutableGraph H = G.copy();
			int[] u_neighbors = H.successorArray(u);
			int u_deg = H.outdegree(u);
			int t_cnt = IntStream.range(0, u_deg).map(i -> {
				int v = u_neighbors[i];
				if (u < v) { // to avoid double counting - for non-optimezed graph

					int[] v_neighbors = H.successorArray(v);
					int v_deg = H.outdegree(v);

					List<Integer> intersection = intersection(u, v, u_neighbors, v_neighbors, u_deg, v_deg);
					int index = intersection.size();

					int degSum = IntStream.range(0, index).map(j -> {
						int w = intersection.get(j); // triangle (u,v,w), where u<v<w

						int[] w_neighbors = H.successorArray(w);
						int w_deg = H.outdegree(w);

						// finds id of triangle
						int id_uvw = findTrinagle_Id(u, v, w);
						int deg4k;

						deg4k = intersect_3(id_uvw,u, v, w, u_neighbors, v_neighbors, w_neighbors, u_deg, v_deg, w_deg);

						_4kDegree[id_uvw] = deg4k;

						return deg4k;
					}).reduce(0, (a, b) -> a + b);
					return index; // this is just the count of triangles
				} else {
					return 0;
				}
			}).reduce(0, (a, b) -> a + b);
			return t_cnt;
		}).mapToObj(BigInteger::valueOf).reduce(BigInteger.ZERO, BigInteger::add);

		long endTime = System.currentTimeMillis();
		System.out.println("Time elapsed (sec) for computing 4kDegrees = " + (endTime - startTime) / 1000.0);

		//String filename = basename + "-" + eta + "-" + "initialsupp_DP_version1.txt";
		// String filename = basename+"eta-"+eta+"initial_Bino.txt";
		// String filename = basename+"eta-"+eta+"initial_poiss_tran.txt";
		//writeResults(_4kDegree, filename);

		System.out.println("total number of triangles: " + t_total_cnt);
		System.out.println("maximum initial 4kdegree: " + max_4kDegree(_4kDegree));
		System.out.println("average initial 4kdegree: " + ave_4kDegree(_4kDegree));

		//ave_4kDegree_nozerotrian(_4kDegree);

	}

	// for each triangle (u,v,w) finds 4-cliques which contains it
	// then, computes probabilistic support using DP
	int intersect_3(int id_uvw, int u, int v, int w, int[] u_neighbors, int[] v_neighbors, int[] w_neighbors, int u_deg, int v_deg,
			int w_deg) {

		int degree = 0; // stores deterministic number of 4-cliques containing a triangle
		int non1prob_count = 0;

		int maxlen = Math.min(u_deg, Math.min(v_deg, w_deg));
		List<Double> _4c_list_prob = new ArrayList<Double>(maxlen);

		for (int i = 0, j = 0, k = 0; i < u_deg && j < v_deg && k < w_deg;) {

			if (u_neighbors[i] == v_neighbors[j] && u_neighbors[i] == w_neighbors[k]) { // Find a 4-clique !
				int z = u_neighbors[i]; // a 4-clique is found

				// check all the triangles have existence prob > eta
				int id_uvz = complex_triangle_Id(u, v, z);
				int id_uwz = complex_triangle_Id(u, w, z);
				int id_vwz = complex_triangle_Id(v, w, z);

				if ((id_uvz != -1) && (id_uwz != -1) && (id_vwz != -1)) { // they have high prob

					int id_uz = Arrays.binarySearch(edgeTail, offset[Math.min(u, z)], offset[Math.min(u, z) + 1],
							Math.max(u, z));
					int id_vz = Arrays.binarySearch(edgeTail, offset[Math.min(v, z)], offset[Math.min(v, z) + 1],
							Math.max(v, z));
					int id_wz = Arrays.binarySearch(edgeTail, offset[Math.min(w, z)], offset[Math.min(w, z) + 1],
							Math.max(w, z));

					double p_uz = edge_prob[id_uz];
					double p_vz = edge_prob[id_vz];
					double p_wz = edge_prob[id_wz];

					double p_exis = p_uz * p_vz * p_wz;

					if (p_exis < epsilon1) { // p_exis != 1
						_4c_list_prob.add(p_exis);
						non1prob_count++;
					}

					degree++;

				}

				i++;
				j++;
				k++;
				continue;
			}

			if (u_neighbors[i] < w_neighbors[k]) {
				if (u_neighbors[i] == v_neighbors[j]) {
					i++;
					j++;
					continue;
				} else if (u_neighbors[i] < v_neighbors[j]) {
					i++;
					continue;
				} else { // u_neighbors[i] > v_neighbors[j]
					j++;
					continue;
				}
			}

			if (u_neighbors[i] > w_neighbors[k]) {
				if (w_neighbors[k] == v_neighbors[j]) {
					k++;
					j++;
					continue;
				} else if (w_neighbors[k] < v_neighbors[j]) {
					k++;
					continue;
				} else { // w_neighbors[k] > v_neighbors[j]
					j++;
					continue;
				}
			}

			if (u_neighbors[i] == w_neighbors[k]) {
				if (u_neighbors[i] < v_neighbors[j]) {
					i++;
					k++;
					continue;
				} else { // v_neighbors[j] is the min
					j++;
					continue;
				}
			}

		}

		int _4kdeg = DP(id_uvw, non1prob_count, _4c_list_prob, degree);

		return _4kdeg;
	}

	// writes the scores in a file
	public void writeResults(int[] support, String filename) throws IOException {
		BufferedWriter wr = new BufferedWriter(new FileWriter(filename));
		for (int i = 0; i < support.length; i++) {
			int u = first_vert[i];
			int v = second_vert[i];
			int w = third_vert[i];
			wr.write(u + "\t" + v + "\t" + w + "\t" + support[i]);
			wr.write("\n");
		}
		wr.close();
	}

	public int find4clique_DPcompute_in_mainPart(int id) {

		int u = first_vert[id];
		int v = second_vert[id];
		int w = third_vert[id];

		int[] u_neigh = G.successorArray(u);
		int[] v_neigh = G.successorArray(v);
		int[] w_neigh = G.successorArray(w);

		int u_deg = G.outdegree(u);
		int v_deg = G.outdegree(v);
		int w_deg = G.outdegree(w);

		// stores the number of 4-cliques, which are alive and contains a triangle
		// the value is without edge probabilities
		int support = 0;
		int non1prob_count = 0;

		int maxlen = Math.min(u_deg, Math.min(v_deg, w_deg));
		List<Double> _4c_list_prob = new ArrayList<Double>(maxlen);

		for (int i = 0, j = 0, k = 0; i < u_deg && j < v_deg && k < w_deg;) {

			if (u_neigh[i] == v_neigh[j] && u_neigh[i] == w_neigh[k]) { // Find a 4-clique !
				int z = u_neigh[i]; // a 4-clique is found

				// now we should check whether triangles (u,v,z), (u,w,z), and (v,w,z)
				// has been processed or not

				int id_uvz = complex_triangle_Id(u, v, z);
				int id_uwz = complex_triangle_Id(u, w, z);
				int id_vwz = complex_triangle_Id(v, w, z);

				//check if they exists
				if ((id_uvz != -1) && (id_uwz != -1) && (id_vwz != -1)) {

					if ((processed.get(id_uvz) == false) && (processed.get(id_uwz) == false)
							&& (processed.get(id_vwz) == false)) {
						int id_uz = Arrays.binarySearch(edgeTail, offset[Math.min(u, z)], offset[Math.min(u, z) + 1],
								Math.max(u, z));
						int id_vz = Arrays.binarySearch(edgeTail, offset[Math.min(v, z)], offset[Math.min(v, z) + 1],
								Math.max(v, z));
						int id_wz = Arrays.binarySearch(edgeTail, offset[Math.min(w, z)], offset[Math.min(w, z) + 1],
								Math.max(w, z));

						double p_uz = edge_prob[id_uz];
						double p_vz = edge_prob[id_vz];
						double p_wz = edge_prob[id_wz];
						double p_exis = p_uz * p_vz * p_wz;

						if (p_exis < epsilon1) { // p_exis != 1

							_4c_list_prob.add(p_exis);
							non1prob_count++;
						}

						support++;

					}
				}
				i++;
				j++;
				k++;
				continue;
			}

			if (u_neigh[i] < w_neigh[k]) {
				if (u_neigh[i] == v_neigh[j]) {
					i++;
					j++;
					continue;
				} else if (u_neigh[i] < v_neigh[j]) {
					i++;
					continue;
				} else { // u_neighbors[i] > v_neighbors[j]
					j++;
					continue;
				}
			}

			if (u_neigh[i] > w_neigh[k]) {
				if (w_neigh[k] == v_neigh[j]) {
					k++;
					j++;
					continue;
				} else if (w_neigh[k] < v_neigh[j]) {
					k++;
					continue;
				} else { // w_neighbors[k] > v_neighbors[j]
					j++;
					continue;
				}
			}

			if (u_neigh[i] == w_neigh[k]) {
				if (u_neigh[i] < v_neigh[j]) {
					i++;
					k++;
					continue;
				} else { // v_neighbors[j] is the min
					j++;
					continue;
				}
			}

		}

		int _4kdeg = DP(id,non1prob_count, _4c_list_prob, support);
		return _4kdeg;

	}

	public int DP(int id_uvw, int n, List<Double> cliques_Exiprobs, int num_of_containingCli) {

		double prob_greater_than_t = 1.0;
		double[][] dp = new double[n + 1][2];

		int sw = 0;
		int eta_support = 0;

		int t = 0;
		while (t <= n && (triangle_prob[id_uvw]*prob_greater_than_t > eta)) {
		//while (t <= n && (prob_greater_than_t > eta)) {

			for (int h = 0; h <= n; h++) { // row
				if (t == 0 && h == 0) {
					dp[h][sw] = 1;
				} else if (h < t) {
					dp[h][sw] = 0;
				} else {

					double prob = cliques_Exiprobs.get(h - 1);
					dp[h][sw] = prob * (t == 0 ? 0 : dp[h - 1][1 - sw]) + (1 - prob) * dp[h - 1][sw];
				}
			}

			double probability = dp[n][sw];

			prob_greater_than_t -= probability; // Pr[sup(e)>=t+1]

			t++;
			sw = 1 - sw;
		}

		if (t == 0) {
			eta_support = num_of_containingCli - n;
		} else {
			eta_support = (num_of_containingCli - n) + t - 1;
		}

		return eta_support;

	}

	// given u,v, and w, finds the id of triangle (u,v,w), where u<v<w
	public int findTrinagle_Id(int target_u, int target_v, int target_w) { // u<v<w

		int edge_uv_ID = Arrays.binarySearch(edgeTail, offset[Math.min(target_u, target_v)],
				offset[Math.min(target_u, target_v) + 1], Math.max(target_u, target_v));

		// now we should look at offset_uv array to see which
		int start_index_search = offset_uv[edge_uv_ID];
		int end_index_search = offset_uv[edge_uv_ID + 1]; // up to this index

		for (int i = start_index_search; i < end_index_search; i++) {
			if (third_vert[i] == target_w) {
				return i;
			}
		}

		return -1;
	}

	List<Integer> intersection(int u, int v, int[] u_neighbors, int[] v_neighbors, int u_deg, int v_deg) {

		int index = 0;
		List<Integer> wlist = new ArrayList<Integer>();
		for (int i = 0, j = 0; i < u_deg && j < v_deg;) {

			if (u_neighbors[i] == v_neighbors[j]) { // Find a triangle !
				// the if condition below is to avoid double counting:
				if (u_neighbors[i] > v) {

					int w = u_neighbors[i];

					int id_uv = Arrays.binarySearch(edgeTail, offset[Math.min(u, v)], offset[Math.min(u, v) + 1],
							Math.max(u, v));
					int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u, w)], offset[Math.min(u, w) + 1],
							Math.max(u, w));
					int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(w, v)], offset[Math.min(w, v) + 1],
							Math.max(w, v));
					double exiprob = edge_prob[id_uv] * edge_prob[id_uw] * edge_prob[id_vw];

					if (exiprob > eta) {
						wlist.add(u_neighbors[i]);
						index++;
					}
				}
				i++;
				j++;
				continue;
			}

			if (u_neighbors[i] < v_neighbors[j]) {
				i++;
				continue;
			}

			if (u_neighbors[i] > v_neighbors[j]) {
				j++;
				continue;
			}
		}

		return wlist;
	}

	int intersection2(int u, int v, int[] u_neighbors, int[] v_neighbors, int u_deg, int v_deg,
			int count_triangles_ddetected_so_far) {

		int index = 0;

		for (int i = 0, j = 0; i < u_deg && j < v_deg;) {

			if (u_neighbors[i] == v_neighbors[j]) { // Find a triangle !

				if (u_neighbors[i] > v) { // to avoid double counting:

					int w = u_neighbors[i]; // triangle (u,v,w)

					int id_uv = Arrays.binarySearch(edgeTail, offset[Math.min(u, v)], offset[Math.min(u, v) + 1],
							Math.max(u, v));
					int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u, w)], offset[Math.min(u, w) + 1],
							Math.max(u, w));
					int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(w, v)], offset[Math.min(w, v) + 1],
							Math.max(w, v));

					double exiprob = edge_prob[id_uv] * edge_prob[id_uw] * edge_prob[id_vw];

					if (exiprob > eta) {
						first_vert[count_triangles_ddetected_so_far + index] = u;
						second_vert[count_triangles_ddetected_so_far + index] = v;
						third_vert[count_triangles_ddetected_so_far + index] = w;
						triangle_prob[count_triangles_ddetected_so_far + index] = exiprob;
						index++;
					}

				}
				i++;
				j++;
				continue;
			}

			if (u_neighbors[i] < v_neighbors[j]) {
				i++;
				continue;
			}

			if (u_neighbors[i] > v_neighbors[j]) {
				j++;
				continue;
			}
		}
		return index;
	}

	private int findHead(int target) {
		// find index such that offset[index] <= target
		// and target < offset[index+1]
		int start = 0;
		int end = offset.length - 1;
		while (start + 1 < end) {
			int mid = start + (end - start) / 2;
			if (offset[mid] <= target) {
				start = mid;
			} else {
				end = mid;
			}
		}
		if (offset[end] <= target) {
			return end;
		}
		if (offset[start] <= target) {
			return start;
		}
		return 0;
	}

	public int max_4kDegree(int[] degree) {
		int max = 0;
		for (int i = 0; i < _4kDegree.length; i++) {
			int val = _4kDegree[i];
			if (val >= max) {
				max = val;
			}
		}
		return max;
	}

	public double ave_4kDegree(int[] degree) {
		double sum = 0.0;
		for (int i = 0; i < degree.length; i++) {
			sum += degree[i];
		}
		double ave = (double) sum / degree.length;
		return ave;
	}

	public void ave_4kDegree_nozerotrian(int[] degree) {
		double sum1 = 0.0;
		int count1 = 0;

		double sum2 = 0.0;
		int count2 = 0;

		for (int i = 0; i < degree.length; i++) {

			// if (!(Exiprob[i] < eta)) {
			// sum1 += degree[i];
			// count1++;
			// }
			if (degree[i] != 0) {
				sum2 += degree[i];
				count2++;
			}

		}
		double ave1 = (double) sum1 / count1;
		double ave2 = (double) sum2 / count2;

		System.out.println("The average score over all the triangles with high existence prob is: " + ave1);
		System.out.println("number of triangles in this case is: " + count1);

		System.out.println("The average score over all the triangles with non-zero final score is: " + ave2);
		System.out.println("number of triangles in this case is: " + count2);

	}

	public void nucleusDecomposition() throws IOException {

		// step 1:sorting triangles based on their support
		System.out.println();
		System.out.println("Starting step 2: sorting the trinagles in ascending order of their 4kdegree.");
		long startTime2 = System.currentTimeMillis();

		int max4kdeg = max_4kDegree(_4kDegree);
		int[] bin = new int[max4kdeg + 1];

		int triangle_count = num_tria_alive;
		for (int i = 0; i < triangle_count; i++) {
			bin[_4kDegree[i]]++;
		}

		int start = 0;
		for (int d = 0; d <= max4kdeg; d++) {
			int num = bin[d];
			bin[d] = start;
			start += num;
		}

		for (int t = 0; t < triangle_count; t++) {
			TrianglePos[t] = bin[_4kDegree[t]];
			sortedTrinagle[TrianglePos[t]] = t;

			bin[_4kDegree[t]]++;
		}

		for (int d = max4kdeg; d >= 1; d--)
			bin[d] = bin[d - 1];
		bin[0] = 0;
		long endTime2 = System.currentTimeMillis();
		System.out.println("Time elapsed (sec) for triangles sorting = " + (endTime2 - startTime2) / 1000.0);

		// main step
		System.out.println("main step of the algorithm");
		long startTime3 = System.currentTimeMillis();

		for (int i = 0; i < triangle_count;) {
			//if (i%1000==0) {
				//System.out.println(i);
			//}

			int id_uvw = sortedTrinagle[i]; // triangle with smallest 4kdegree
			if (processed.get(id_uvw) == false) {
				int current_deg = _4kDegree[id_uvw];

				// we assume that u,v,w have been already stored. So, using id of the triangle
				// we can find u,v,w, where u<v<w.
				int u = first_vert[id_uvw];
				int v = second_vert[id_uvw];
				int w = third_vert[id_uvw];

				// 4clique finding
				int[] u_neighbors = G.successorArray(u);
				int u_deg = G.outdegree(u);
				int[] v_neighbors = G.successorArray(v);
				int v_deg = G.outdegree(v);
				int[] w_neighbors = G.successorArray(w);
				int w_deg = G.outdegree(w);

				List<Integer> clique_list = Find_4clique(u, v, w, u_neighbors, v_neighbors, w_neighbors, u_deg, v_deg,
						w_deg);
				int _4kdeg = clique_list.size();

				processed.set(id_uvw);

				for (int j = 0; j < _4kdeg; j++) {
					int z = clique_list.get(j); // 4clique: (u,v,w,z)

					// all triangles in the 4clique
					int id_uvz = complex_triangle_Id(u, v, z);
					int id_uwz = complex_triangle_Id(u, w, z);
					int id_vwz = complex_triangle_Id(v, w, z);

					if ((id_uvz != -1) && (id_uwz != -1) && (id_vwz != -1)) {

						if ((processed.get(id_uvz) == false) && (processed.get(id_uwz) == false)
								&& (processed.get(id_vwz) == false)) {

							int curDeg_id_uvz = _4kDegree[id_uvz];
							if (curDeg_id_uvz > current_deg) {

								int updatedDeg_id_uvz = find4clique_DPcompute_in_mainPart(id_uvz);

								if (Math.abs(updatedDeg_id_uvz - curDeg_id_uvz) > 1) {
									updatedDeg_id_uvz = Math.max(curDeg_id_uvz - 1, 0);
								} // if numerical instability happens we fix it

								if (Math.abs(updatedDeg_id_uvz - curDeg_id_uvz) == 1) {
									// decrement and swap it in the sortedTrinagle array
									int pos_id_uvz = TrianglePos[id_uvz];
									int pos_tri_same_bin = bin[curDeg_id_uvz];
									int id_tri_same_bin = sortedTrinagle[pos_tri_same_bin];

									if (id_uvz != id_tri_same_bin) {
										TrianglePos[id_uvz] = pos_tri_same_bin;
										sortedTrinagle[pos_id_uvz] = id_tri_same_bin;
										TrianglePos[id_tri_same_bin] = pos_id_uvz;
										sortedTrinagle[pos_tri_same_bin] = id_uvz;
									}
									bin[curDeg_id_uvz]++;
									_4kDegree[id_uvz]--;
								}

							}

							int curDeg_id_uwz = _4kDegree[id_uwz];
							if (curDeg_id_uwz > current_deg) {
								int updatedDeg_id_uwz = find4clique_DPcompute_in_mainPart(id_uwz);

								if (Math.abs(updatedDeg_id_uwz - curDeg_id_uwz) > 1) {
									updatedDeg_id_uwz = Math.max(curDeg_id_uwz - 1, 0);
								} // we fix it, if numerical stability happens

								if (Math.abs(updatedDeg_id_uwz - curDeg_id_uwz) == 1) {
									int pos_id_uwz = TrianglePos[id_uwz];
									int pos_tri_same_bin = bin[curDeg_id_uwz];
									int id_tri_same_bin = sortedTrinagle[pos_tri_same_bin];

									if (id_uwz != id_tri_same_bin) {
										TrianglePos[id_uwz] = pos_tri_same_bin;
										sortedTrinagle[pos_id_uwz] = id_tri_same_bin;
										TrianglePos[id_tri_same_bin] = pos_id_uwz;
										sortedTrinagle[pos_tri_same_bin] = id_uwz;
									}
									bin[curDeg_id_uwz]++;
									_4kDegree[id_uwz]--;
								}

							}

							int curDeg_id_vwz = _4kDegree[id_vwz];
							if (curDeg_id_vwz > current_deg) {
								int updatedDeg_id_vwz = find4clique_DPcompute_in_mainPart(id_vwz);

								if (Math.abs(updatedDeg_id_vwz - curDeg_id_vwz) > 1) {
									// numerical instability
									updatedDeg_id_vwz = Math.max(curDeg_id_vwz - 1, 0);
								}

								if (Math.abs(updatedDeg_id_vwz - curDeg_id_vwz) == 1) {
									int pos_id_vwz = TrianglePos[id_vwz];
									int pos_tri_same_bin = bin[curDeg_id_vwz];
									int id_tri_same_bin = sortedTrinagle[pos_tri_same_bin];

									if (id_vwz != id_tri_same_bin) {
										TrianglePos[id_vwz] = pos_tri_same_bin;
										sortedTrinagle[pos_id_vwz] = id_tri_same_bin;
										TrianglePos[id_tri_same_bin] = pos_id_vwz;
										sortedTrinagle[pos_tri_same_bin] = id_vwz;
									}
									bin[curDeg_id_vwz]++;
									_4kDegree[id_vwz]--;
								}

							} // if for id_vwz

						}
					}

				}

			}
			i++;
		}

		long endTime3 = System.currentTimeMillis();
		System.out.println("Time elapsed (sec) for main part= " + (endTime3 - startTime3) / 1000.0);

		System.out.println("maximum final support score: " + max_4kDegree(_4kDegree));
		System.out.println("average final support score: " + ave_4kDegree(_4kDegree));

		//ave_4kDegree_nozerotrian(_4kDegree);

		String filename = basename + "-" + eta + "-" + "finalsupp_DP_Modified.txt";
		writeResults(_4kDegree, filename);
		
		//for (int i=0;i<_4kDegree.length;i++) {
			//System.out.println(first_vert[i] + " " + second_vert[i] + " " + third_vert[i]+ " " + _4kDegree[i]);
		//}
	}

	List<Integer> Find_4clique(int u, int v, int w, int[] u_neighbors, int[] v_neighbors, int[] w_neighbors, int u_deg,
			int v_deg, int w_deg) {

		int degree = 0;

		int minlength = 0;
		minlength = Math.min(u_deg, Math.min(v_deg, w_deg));
		List<Integer> clique_list = new ArrayList<Integer>(minlength);

		for (int i = 0, j = 0, k = 0; i < u_deg && j < v_deg && k < w_deg;) {

			if (u_neighbors[i] == v_neighbors[j] && u_neighbors[i] == w_neighbors[k]) { // Find a 4-clique !
//	            System.out.println(u + ", " + v + ", " + w  + ", " + u_neighbors[i]); // for checking
				clique_list.add(u_neighbors[i]);
				degree++;
				i++;
				j++;
				k++;
				continue;
			}

			if (u_neighbors[i] < w_neighbors[k]) {
				if (u_neighbors[i] == v_neighbors[j]) {
					i++;
					j++;
					continue;
				} else if (u_neighbors[i] < v_neighbors[j]) {
					i++;
					continue;
				} else { // u_neighbors[i] > v_neighbors[j]
					j++;
					continue;
				}
			}

			if (u_neighbors[i] > w_neighbors[k]) {
				if (w_neighbors[k] == v_neighbors[j]) {
					k++;
					j++;
					continue;
				} else if (w_neighbors[k] < v_neighbors[j]) {
					k++;
					continue;
				} else { // w_neighbors[k] > v_neighbors[j]
					j++;
					continue;
				}
			}

			if (u_neighbors[i] == w_neighbors[k]) {
				if (u_neighbors[i] < v_neighbors[j]) {
					i++;
					k++;
					continue;
				} else { // v_neighbors[j] is the min
					j++;
					continue;
				}
			}

		}
		return clique_list;
	}

	public int complex_triangle_Id(int target_u, int target_v, int target_w) {
		int u;
		int v;
		int w;

		if ((target_u < target_v) && (target_u < target_w)) {
			u = target_u;
			if (target_v < target_w) {
				v = target_v;
				w = target_w;
			} else {
				v = target_w;
				w = target_v;
			}

			// List<Integer> l1 = new ArrayList<>(3);
			// l1.add(u);
			// l1.add(v);
			//// l1.add(w);
			// return bm.getKey(l1);
			return findTrinagle_Id(u, v, w);
		} else if ((target_v < target_u) && (target_v < target_w)) {
			u = target_v;
			if (target_u < target_w) {
				v = target_u;
				w = target_w;
			} else {
				v = target_w;
				w = target_u;
			}
			// List<Integer> l2 = new ArrayList<>(3);
			// l2.add(u);
			// l2.add(v);
			// l2.add(w);
			// return bm.getKey(l2);
			return findTrinagle_Id(u, v, w);
		} else {
			u = target_w;
			if (target_u < target_v) {
				v = target_u;
				w = target_v;
			} else {
				v = target_v;
				w = target_u;
			}
			// List<Integer> l3 = new ArrayList<>(3);
			// l3.add(u);
			// l3.add(v);
			// l3.add(w);
			// return bm.getKey(l3);
			return findTrinagle_Id(u, v, w);
		}

	}

	public static void main(String[] args) throws Exception {
		//long beforeUsedMem=Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory(); 
		
		
		String basename = args[0];
		int precision = Integer.valueOf(args[1]);
		double eta = Double.valueOf(args[2]);
		long startTime = System.currentTimeMillis();
		//long memory = 0;
		//long memory2 = 0;
		//int[][] dp = new int[3][2];
		//System.out.println("length" + dp.length);
		//AtomicIntegerArray scheduled = new AtomicIntegerArray(10); 
		//int[] arr = new int[10];
		//int[] arr2 = new int[20];
		//memory += 10*Integer.BYTES;
		//memory += 20*Integer.SIZE;
		//System.out.println(memory);
		
		//String basename = "biomine-proc.w";
		//String basename = "DBLP_comput-proc.w";
		//int precision = 16;
		//double eta = 0.00000000001;
		//double eta = 0.0000000000001;
		//double eta = 0.1;
		//double eta = 0.000001;
		//double eta = 0.0000000000000001;
		//for (int i=1;i<=1000;i++) {
			//String basename = "test_global_1000"+i+".w";
		    //String basename = "DBLP_data_net_proc.w";
		//String basename = "DBLP_data_net-proc.w";
			PNucl_DP_Efficient pr_nc = new PNucl_DP_Efficient(basename, precision, eta);
			pr_nc.computeTriangles_degree();
			 pr_nc.nucleusDecomposition();
			 String filename = basename + "-" + eta + "-" + "finalsupp_DP_Modified.txt";
			 pr_nc.writeResults(pr_nc._4kDegree, filename);
			 
			 long endTime = System.currentTimeMillis();
			 System.out.println("Total running time is = " + (endTime - startTime) / 1000.0);
		//}
		//String basename = "test_global_10001.w";
		// String basename = "Flickr-proc.w";
		//String basename = "testgraph1.w";
		//long tircount = 7;
		// long tircount = 4582169;
		

		// double eta = 0.1;
		// String basename = "newTest.w";
		// long tircount = 6;
		// int precision = 2;
		// double eta = 0.0;

		

		//PNucl_DP_Efficient pr_nc = new PNucl_DP_Efficient(basename, precision, eta);

		/*for (int i = 0; i < pr_nc_DPv1.first_vert.length; i++) {
			System.out.println(
					pr_nc_DPv1.first_vert[i] + " " + pr_nc_DPv1.second_vert[i] + " " + pr_nc_DPv1.third_vert[i]);
		}
		System.out.println(pr_nc_DPv1.complex_triangle_Id(3, 4, 1));*/
		// System.out.println(pr_nc_DPv1.findTrinagle_Id(0, 2, 3));
		// System.out.println(pr_nc_DPv1.findTrinagle_Id(1, 2, 3));
		// System.out.println(pr_nc_DPv1.findTrinagle_Id(2, 3, 4));

		//pr_nc.computeTriangles_degree();
		// pr_nc.nucleusDecomposition();
		 
		// long afterUsedMem=Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
		// long actualMemUsed=afterUsedMem-beforeUsedMem;
		// System.out.println("the amount of memory used " + actualMemUsed);
		// for (int i=0;i<pr_nc_DPv1._4kDegree.length;i++) {
		// System.out.println("trinagle is:" + pr_nc_DPv1.first_vert[i] + " " +
		// pr_nc_DPv1.second_vert[i] + " " + pr_nc_DPv1.third_vert[i]);
		// System.out.println(pr_nc_DPv1._4kDegree[i]);
		// }

	}
}
