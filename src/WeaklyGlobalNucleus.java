
import it.unimi.dsi.webgraph.ImmutableGraph;
import it.unimi.dsi.webgraph.labelling.ArcLabelledImmutableGraph;
import it.unimi.dsi.webgraph.labelling.Label;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import java.util.Set;
import java.util.SortedMap;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;

import java.util.stream.IntStream;

import gnu.trove.map.TIntObjectMap;

import gnu.trove.map.hash.TIntObjectHashMap;

public class WeaklyGlobalNucleus {

	String basename;
	int precision;
	double eta;

	int num_tria_alive;
	int max_score;

	BitSet visited;
	int[] support; // stores 4kdegree for each triangle

	String filename;
	int num_samples;

	ArcLabelledImmutableGraph G; // uncertain graph
	int n; // vertices of G
	int m; // edges of G

	final int[] edgeTail;
	final int[] offset; // offset array for edges
	double[] edge_prob;

	BitSet member_of_component;
	BitSet member_of_component_for_edge;

	// for each triangle (u,v,w), where u<v<w:
	int[] first_vert; // first_vert stores u
	int[] second_vert; // second_vert stores v
	int[] third_vert; // third_vert stores
	int[] offset_uv; // offset array for triangles

	Map<Integer, List<Integer>> triangle_comneigh;
	Map<Integer, Integer> Degree_in_sample;

	int[] glob_score;
	BitSet edge_in_sample; // it says whether an edge appear in the sample or not

	BigInteger t_total_cnt = BigInteger.ZERO;
	BitSet processed; // records if a triangle has been processed or not

	List<Integer> list_of_triangles;
	List<Integer> list_of_vertices;
	List<Integer> list_of_edges;

	public WeaklyGlobalNucleus(String basename, int precision, double eta, int num_samples)
			throws Exception {

		this.basename = basename;
		this.precision = precision;
		this.eta = eta;
		this.num_samples = num_samples;

		this.filename = this.basename + "-" + this.eta + "-" + "finalsupp_DP_Modified.txt";

		System.out.println("loading the graph...");
		long startTime1 = System.currentTimeMillis();
		G = ArcLabelledImmutableGraph.load(basename);
		n = G.numNodes();
		m = (int) (G.numArcs() / 2);
		long endTime1 = System.currentTimeMillis();
		System.out.println("Time elapsed (sec) for loading the graph = " + (endTime1 - startTime1) / 1000.0);

		edgeTail = new int[m]; // store tail for each edge
		offset = new int[n + 1]; // edge offset for edgeHead array
		edge_prob = new double[m];

		member_of_component = new BitSet(n);
		member_of_component_for_edge = new BitSet(m);

		System.out.println("Start storing triangles in edgeTail...");
		long startTime2 = System.currentTimeMillis();
		Populate_edges();
		long endTime2 = System.currentTimeMillis();
		System.out.println("Time elapsed (sec) for storing edges = " + (endTime2 - startTime2) / 1000.0);

		System.out.println("start enumerating triangles...");
		long startTime3 = System.currentTimeMillis();
		EnumerateTriangles(); // set the value of "num_tria_alive"
		long endTime3 = System.currentTimeMillis();
		System.out.println("Time elapsed (sec) for enumeration = " + (endTime3 - startTime3) / 1000.0);

		System.out.println("the number of triangles with high existence probability is: " + this.num_tria_alive);

		first_vert = new int[num_tria_alive];
		second_vert = new int[num_tria_alive];
		third_vert = new int[num_tria_alive];
		offset_uv = new int[m + 1]; // starting index of triangles with edge (u,v) where u<v
		support = new int[num_tria_alive]; // store nucleus score for each triangle
		visited = new BitSet(num_tria_alive);

		long startTime4 = System.currentTimeMillis();
		System.out.println("start storing triangles...");
		Populate_triangles();
		long endTime4 = System.currentTimeMillis();
		System.out.println("Time elapsed (sec) for storing trinagles = " + (endTime4 - startTime4) / 1000.0);

		long startTime5 = System.currentTimeMillis();
		System.out.println("start reading nucleus scores...");
		Readscore();
		long endTime5 = System.currentTimeMillis();
		System.out.println("Time elapsed (sec) for reading nucleus scores = " + (endTime5 - startTime5) / 1000.0);

		int max = 0;
		for (int i = 0; i < support.length; i++) {
			int val = support[i];
			if (val >= max) {
				max = val;
			}
		}

		max_score = max;
		System.out.println("maximum nucleuss score is: " + max_score);

		list_of_triangles = new ArrayList<Integer>(num_tria_alive);
		list_of_vertices = new ArrayList<Integer>(n);
		list_of_edges = new ArrayList<Integer>(m);

		triangle_comneigh = new HashMap<>();
		Degree_in_sample = new HashMap<>();

		glob_score = new int[num_tria_alive];
		edge_in_sample = new BitSet(m);

		processed = new BitSet(num_tria_alive); // no triangle has been processed yet

	}
	public double find_clust_coef(List<Integer> vert_list, List<Integer> edge_list, List<Integer> tri_list) {

		double sum_numer = 0.0;
		double sum_denom = 0.0;
		
		BitSet vert = new BitSet(n);
		BitSet edge = new BitSet(m);
		
		for (int i=0;i<vert_list.size();i++) {
			int u = vert_list.get(i);
			vert.set(u);
		}
		
		for (int i=0;i<edge_list.size();i++) {
			int e = edge_list.get(i);
			edge.set(e);
		}
		
		//System.out.println(vert.cardinality());
		//System.out.println(edge.cardinality() + " " + edge_list.size());
		
		
		
		

		for (int i=0;i<vert_list.size();i++) {
			int u = vert_list.get(i);
			int[] u_neigh = G.successorArray(u);
			int len_u_neigh = u_neigh.length;
			for (int j = 0; j < len_u_neigh; j++) {
				int neigh_1 = u_neigh[j];
				if (vert.get(neigh_1) == true) {
					int id_1 = Arrays.binarySearch(edgeTail, offset[Math.min(u, neigh_1)],
							offset[Math.min(u, neigh_1) + 1], Math.max(u, neigh_1));
					if (edge.get(id_1)==true) { // (u,neigh1) exists in the graph
						for (int k = j + 1; k < len_u_neigh; k++) {
							int neigh_2 = u_neigh[k];
							if (vert.get(neigh_2) == true) {
								int id_2 = Arrays.binarySearch(edgeTail, offset[Math.min(u, neigh_2)],
										offset[Math.min(u, neigh_2) + 1], Math.max(u, neigh_2));
								if (edge.get(id_2) == true) {
									double prob = edge_prob[id_1] * edge_prob[id_2];
									sum_denom += prob;
								}
							}
						}

					}
				}
			}
		}

		for (int i = 0; i < tri_list.size(); i++) {
			int id_uvw = tri_list.get(i);
			int u = first_vert[id_uvw];
			int v = second_vert[id_uvw];
			int w = third_vert[id_uvw];

			int id_uv = Arrays.binarySearch(edgeTail, offset[Math.min(u, v)], offset[Math.min(u, v) + 1],
					Math.max(u, v));
			int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u, w)], offset[Math.min(u, w) + 1],
					Math.max(u, w));
			int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(v, w)], offset[Math.min(v, w) + 1],
					Math.max(v, w));

			double Exiprob = edge_prob[id_uv] * edge_prob[id_uw] * edge_prob[id_vw];
			sum_numer += Exiprob;
		}

		double cc = (double) (3 * sum_numer) / (sum_denom);
		return cc;
	}

	public double find_clust_coef(Set<Integer> vert_list, Set<Integer> edge_list, List<Integer> tri_list) {

		double sum_numer = 0.0;
		double sum_denom = 0.0;

		for (Integer u : vert_list) {
			int[] u_neigh = G.successorArray(u);
			int len_u_neigh = u_neigh.length;
			for (int j = 0; j < len_u_neigh; j++) {
				int neigh_1 = u_neigh[j];
				if (vert_list.contains(neigh_1)) {
					int id_1 = Arrays.binarySearch(edgeTail, offset[Math.min(u, neigh_1)],
							offset[Math.min(u, neigh_1) + 1], Math.max(u, neigh_1));
					if (edge_list.contains(id_1)) { // (u,neigh1) exists in the graph
						for (int k = j + 1; k < len_u_neigh; k++) {
							int neigh_2 = u_neigh[k];
							if (vert_list.contains(neigh_2)) {
								int id_2 = Arrays.binarySearch(edgeTail, offset[Math.min(u, neigh_2)],
										offset[Math.min(u, neigh_2) + 1], Math.max(u, neigh_2));
								if (edge_list.contains(id_2)) {
									double prob = edge_prob[id_1] * edge_prob[id_2];
									sum_denom += prob;
								}
							}
						}

					}
				}
			}
		}

		for (int i = 0; i < tri_list.size(); i++) {
			int id_uvw = tri_list.get(i);
			int u = first_vert[id_uvw];
			int v = second_vert[id_uvw];
			int w = third_vert[id_uvw];

			int id_uv = Arrays.binarySearch(edgeTail, offset[Math.min(u, v)], offset[Math.min(u, v) + 1],
					Math.max(u, v));
			int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u, w)], offset[Math.min(u, w) + 1],
					Math.max(u, w));
			int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(v, w)], offset[Math.min(v, w) + 1],
					Math.max(v, w));

			double Exiprob = edge_prob[id_uv] * edge_prob[id_uw] * edge_prob[id_vw];
			sum_numer += Exiprob;
		}

		double cc = (double) (3 * sum_numer) / (sum_denom);
		return cc;
	}

	public double find_density(List<Integer> edge_list, List<Integer> vert_list) {
		double sum = 0.0;
		for (int e = 0; e < edge_list.size(); e++) {
			int id_e = edge_list.get(e);
			sum += edge_prob[id_e];
		}

		int num_vertices = vert_list.size();
		return sum / ((double) 0.5 * num_vertices * (num_vertices - 1));
	}

	public Map<String, Integer> Read_author_id(String filename1) throws IOException {
		Map<String, Integer> author_name_ID = new HashMap<>();
		BufferedReader br2 = new BufferedReader(new FileReader(filename1));
		String line2 = br2.readLine();

		int count = 0;
		while (line2 != null) {
			int id = Integer.parseInt(line2.split("\\s+")[0]);
			String author = line2.split("\\s+")[1];
			author_name_ID.put(author, id);
			count++;

			line2 = br2.readLine();
		}

		return author_name_ID;

	}

	// outputs the number of triangles greater than eta (threshold)
	public void EnumerateTriangles() throws Exception {

		t_total_cnt = IntStream.range(0, n).map(u -> {
			// if (u % 1_000_000 == 0)
			// System.out.println(u);

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

		// System.out.println("Number of Triangles with high existence probability: " +
		// t_total_cnt);
		this.num_tria_alive = t_total_cnt.intValue();
	}

	List<Integer> intersection(int u, int v, int[] u_neighbors, int[] v_neighbors, int u_deg, int v_deg) {

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

	public void Readscore() throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String line = br.readLine();
		int count = 0;

		while (line != null) {

			support[count] = Integer.parseInt(line.split("\t")[3]);
			count++;
			line = br.readLine();
		}
		br.close();
	}

	public void DFS_iter(int id, int k_val) {

		Stack<Integer> stack = new Stack<>();
		stack.push(id);

		while (!(stack.isEmpty())) {
			id = stack.peek();
			stack.pop();

			int u = first_vert[id];
			int v = second_vert[id];
			int w = third_vert[id];

			if (visited.get(id) == false) {

				visited.set(id);
				list_of_triangles.add(id);
				int id_uv = Arrays.binarySearch(edgeTail, offset[Math.min(u, v)], offset[Math.min(u, v) + 1],
						Math.max(u, v));
				int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u, w)], offset[Math.min(u, w) + 1],
						Math.max(u, w));
				int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(v, w)], offset[Math.min(v, w) + 1],
						Math.max(v, w));

				if (member_of_component_for_edge.get(id_uv) == false) {
					list_of_edges.add(id_uv);
					member_of_component_for_edge.set(id_uv);
				}
				if (member_of_component_for_edge.get(id_uw) == false) {
					list_of_edges.add(id_uw);
					member_of_component_for_edge.set(id_uw);
				}
				if (member_of_component_for_edge.get(id_vw) == false) {
					list_of_edges.add(id_vw);
					member_of_component_for_edge.set(id_vw);
				}

				if (member_of_component.get(u) == false) {
					list_of_vertices.add(u);
					member_of_component.set(u);
				}
				if (member_of_component.get(v) == false) {
					list_of_vertices.add(v);
					member_of_component.set(v);
				}
				if (member_of_component.get(w) == false) {

					list_of_vertices.add(w);
					member_of_component.set(w);
				}

				int[] u_neighbors = G.successorArray(u);
				int[] v_neighbors = G.successorArray(v);
				int[] w_neighbors = G.successorArray(w);
				int u_deg = G.outdegree(u);
				int v_deg = G.outdegree(v);
				int w_deg = G.outdegree(w);

				List<Integer> _4cli_list = Find_4clique(u, v, w, u_neighbors, v_neighbors, w_neighbors, u_deg, v_deg,
						w_deg);

				int len = _4cli_list.size();
				List<Integer> neigh = new ArrayList<>(len * 3);

				for (int j = 0; j < len; j++) {
					int z = _4cli_list.get(j);

					int id_uvz = complex_triangle_Id(u, v, z);
					int id_uwz = complex_triangle_Id(u, w, z);
					int id_vwz = complex_triangle_Id(v, w, z);

					if ((id_uvz != -1) && (id_uwz != -1) && (id_vwz != -1)) {
						if ((support[id_uvz] >= k_val) && (support[id_uwz] >= k_val) && (support[id_vwz] >= k_val)) {
							neigh.add(id_uvz);
							neigh.add(id_uwz);
							neigh.add(id_vwz);
						}
					}
				} // we find all the neighbors

				for (int i = 0; i < neigh.size(); i++) {
					int id_neigh = neigh.get(i);
					if (visited.get(id_neigh) == false) {

						stack.push(id_neigh);
					}
				}

			}
		}
	}

	public double findDensity_and_edges_vertices(List<Integer> list_vert, List<Integer> list_edge) {

		int len = list_vert.size(); // number of vertices
		int len_edge = list_edge.size(); // number of edges

		double sum = 0.0;

		for (int i = 0; i < len_edge; i++) {
			int id_edge = list_edge.get(i);
			double p_edge = edge_prob[id_edge];
			sum += p_edge;
		}

		int num_vertices = len;
		// int num_edges = len_edge;

		double desity_value = sum / ((double) 0.5 * num_vertices * (num_vertices - 1));
		// System.out.println("density " + density);

		return desity_value;
	}

	public double find_clust_coeif(List<Integer> list_vert, List<Integer> list_edge) {

		int len = list_vert.size();

		double sum_numer = 0.0;
		double sum_denom = 0.0;

		for (int i = 0; i < len; i++) {
			int u = list_vert.get(i); // vertex in the list
			int[] u_neigh = G.successorArray(u);
			int len_u_neigh = u_neigh.length;

			for (int j = 0; j < len_u_neigh; j++) {
				int neigh_1 = u_neigh[j];
				if (member_of_component.get(neigh_1) == true) {
					for (int k = j + 1; k < len_u_neigh; k++) {
						int neigh_2 = u_neigh[k];

						if ((member_of_component.get(neigh_2) == true)) {

							int id_uneigh1 = Arrays.binarySearch(edgeTail, offset[Math.min(u, neigh_1)],
									offset[Math.min(u, neigh_1) + 1], Math.max(u, neigh_1));
							int id_uneigh2 = Arrays.binarySearch(edgeTail, offset[Math.min(u, neigh_2)],
									offset[Math.min(u, neigh_2) + 1], Math.max(u, neigh_2));

							if ((member_of_component_for_edge.get(id_uneigh1)) == true
									&& (member_of_component_for_edge.get(id_uneigh2) == true)) {
								double prob = edge_prob[id_uneigh1] * edge_prob[id_uneigh2];
								sum_denom += prob;
							}

						}
					}
				}
			}

		}

		// Find triangles

		for (int e = 0; e < m; e++) {

			if (member_of_component_for_edge.get(e) == true) {
				int u = findHead(e); // u<v
				int v = edgeTail[e];

				int[] u_neighbors = G.successorArray(u);
				int[] v_neighbors = G.successorArray(v);

				int u_deg = G.outdegree(u);
				int v_deg = G.outdegree(v);

				for (int i = 0, j = 0; i < u_deg && j < v_deg;) {

					if (u_neighbors[i] == v_neighbors[j]) { // Find a triangle !
						// the if condition below is to avoid double counting:
						if (u_neighbors[i] > v) {

							int w = u_neighbors[i];

							int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u, w)],
									offset[Math.min(u, w) + 1], Math.max(u, w));
							int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(w, v)],
									offset[Math.min(w, v) + 1], Math.max(w, v));

							if ((member_of_component_for_edge.get(id_uw) == true)
									&& (member_of_component_for_edge.get(id_vw) == true)) {
								double Exiprob = edge_prob[e] * edge_prob[id_uw] * edge_prob[id_vw];
								sum_numer += Exiprob;
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

			}
		}

		// System.out.println(sum_denom);
		double cc = (double) (3 * sum_numer) / (sum_denom);
		return cc;

	}

	public void Main_BasicSampling() throws Exception {
		
		

		String wglobalfile = "WeaklyGlobalNucleus-" + num_samples + "-" + basename + "-" + eta + "-.txt";
		BufferedWriter wg = new BufferedWriter(new FileWriter(wglobalfile));

		//String localfile = "Local-" + num_samples + "-" + basename + "-" + eta + "-.txt";
		//BufferedWriter l = new BufferedWriter(new FileWriter(localfile));

		long startTime1 = System.currentTimeMillis();
		System.out.println("start finding connected components and their weakly-global subgraph....");

		for (int k = max_score; k > 0; k--) {
			//if (k==1) {
				//break;
			//}
			int count_con = 0;
			visited.clear();
			wg.write("k-nucleus " + k + "=" + "\n");

			for (int id = 0; id < num_tria_alive; id++) {
				if ((visited.get(id) == false) && (support[id] >= k)) {

					list_of_vertices = new ArrayList<>();
					list_of_triangles = new ArrayList<>();
					list_of_edges = new ArrayList<>();

					member_of_component.clear();
					member_of_component_for_edge.clear();

					DFS_iter(id, k);
					
					//double density_local = find_density(list_of_edges, list_of_vertices);
					//double cc_local = find_clust_coef(list_of_vertices, list_of_edges, list_of_triangles);
					//l.write(density_local + "\t" + cc_local + "\t" +  list_of_edges.size() + "\t" + list_of_vertices.size() + "\n");
					
					// System.out.println(list_of_triangles.size());
					count_con++;

					// Perform weakly-global decomposition

					SortedMap<Integer, SortedMap<Integer, Double>> G_conn = new TreeMap<>();

					for (int i = 0; i < list_of_edges.size(); i++) {
						int id_e = list_of_edges.get(i);
						int u = findHead(id_e);
						int v = edgeTail[id_e];

						if (!G_conn.containsKey(u)) {
							G_conn.put(u, new TreeMap<>());
						}

						if (!G_conn.containsKey(v)) {
							G_conn.put(v, new TreeMap<>());
						}

						G_conn.get(u).put(v, edge_prob[id_e]);
						G_conn.get(v).put(u, edge_prob[id_e]);
					}

					Arrays.fill(glob_score, 0); // before sampling for each component

					for (int N = 0; N < num_samples; N++) {

						// neighbor set for each triangle in the sample
						triangle_comneigh.clear();

						// degree for each triangle in the sample
						Degree_in_sample.clear();

						if (N % 100 == 0) {
							System.out.println("the sample " + N + " for k = " + k);
						}

						edge_in_sample.clear();
						// Arrays.fill(_4kDegree, -1);
						processed.clear();
						// triangle_comneigh.clear();
						// triangle_clique.clear();

						sampling(); // set edge_in_sample for the edges alive in sample.
						// System.out.println("number of edges in the sample: " +

						// System.out.println(member_of_component_for_edge.cardinality() + " " +
						// edge_in_sample.cardinality());

						// now, compute k-nucleus decomposition

						// 1) find triangles alive and their degrees

						for (Integer key_u : G_conn.keySet()) {
							for (Integer key_v : G_conn.get(key_u).keySet()) {
								int u = key_u.intValue();
								int v = key_v.intValue();
								if (u < v) {

									int id_uv = Arrays.binarySearch(edgeTail, offset[Math.min(u, v)],
											offset[Math.min(u, v) + 1], Math.max(u, v));

									if (edge_in_sample.get(id_uv) == true) {

										ArrayList<Integer> u_neigh = new ArrayList<Integer>(G_conn.get(key_u).keySet());
										ArrayList<Integer> v_neigh = new ArrayList<Integer>(G_conn.get(key_v).keySet());
										int u_deg = u_neigh.size();
										int v_deg = v_neigh.size();

										List<Integer> tri_list = intersection_on_SmallGraph(u, v, id_uv, u_neigh,
												v_neigh, u_deg, v_deg);
										int index = tri_list.size();
										for (int j = 0; j < index; j++) {
											int w = tri_list.get(j);
											int id_uvw = findTrinagle_Id(u, v, w);
											ArrayList<Integer> w_neigh = new ArrayList<Integer>(G_conn.get(w).keySet());
											int w_deg = w_neigh.size();
											List<Integer> arr = intersect_3_on_samllGraph(u, v, w, u_neigh, v_neigh,
													w_neigh, u_deg, v_deg, w_deg, id_uvw);

											triangle_comneigh.put(id_uvw, arr);
										}
									}
								}
							}
						}

						int max4kdeg = 0;
						for (Integer Id_uvw : triangle_comneigh.keySet()) {
							int deg = triangle_comneigh.get(Id_uvw).size();
							Degree_in_sample.put(Id_uvw, deg);
							if (deg >= max4kdeg) {
								max4kdeg = deg;
							}
						}

						@SuppressWarnings("unchecked")
						Set<Integer>[] D = new Set[max4kdeg + 1];

						for (Integer id_uvw : Degree_in_sample.keySet()) {
							int deg = Degree_in_sample.get(id_uvw);
							if (D[deg] == null) {
								D[deg] = new TreeSet<Integer>();
							}
							D[deg].add(id_uvw);
						}

						for (int deg = 0; deg < D.length; deg++) {

							Set<Integer> s = D[deg];

							if (deg > k) {
								break;
							}

							while (s != null && !s.isEmpty()) {

								int id_uvw = ((TreeSet<Integer>) s).first();
								s.remove(id_uvw);

								Degree_in_sample.replace(id_uvw, deg);

								int u = first_vert[id_uvw];
								int v = second_vert[id_uvw];
								int w = third_vert[id_uvw];

								List<Integer> clique_list = triangle_comneigh.get(id_uvw);
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

											int PrevDeg_id_uvz = Degree_in_sample.get(id_uvz);
											if (PrevDeg_id_uvz > deg) {
												int newDeg_id_uvz = Math.max(PrevDeg_id_uvz - 1, 0);
												D[PrevDeg_id_uvz].remove(id_uvz);
												if (D[newDeg_id_uvz] == null) {
													D[newDeg_id_uvz] = new TreeSet<Integer>();
												}
												D[newDeg_id_uvz].add(id_uvz);
												Degree_in_sample.replace(id_uvz, newDeg_id_uvz);
											}

											int PrevDeg_id_uwz = Degree_in_sample.get(id_uwz);
											if (PrevDeg_id_uwz > deg) {
												int newDeg_id_uwz = Math.max(PrevDeg_id_uwz - 1, 0);
												D[PrevDeg_id_uwz].remove(id_uwz);
												if (D[newDeg_id_uwz] == null) {
													D[newDeg_id_uwz] = new TreeSet<Integer>();
												}
												D[newDeg_id_uwz].add(id_uwz);
												Degree_in_sample.replace(id_uwz, newDeg_id_uwz);
											}

											int PrevDeg_id_vwz = Degree_in_sample.get(id_vwz);
											if (PrevDeg_id_vwz > deg) {
												int newDeg_id_vwz = Math.max(PrevDeg_id_vwz - 1, 0);
												D[PrevDeg_id_vwz].remove(id_vwz);
												if (D[newDeg_id_vwz] == null) {
													D[newDeg_id_vwz] = new TreeSet<Integer>();
												}
												D[newDeg_id_vwz].add(id_vwz);
												Degree_in_sample.replace(id_vwz, newDeg_id_vwz);
											}
										}
									}
								} // end checking neighbors
							}
						} // done decomposition

						// For each sample, check the triangles which are part of a k-nucleus

						for (Integer Id_uvw : Degree_in_sample.keySet()) {
							int val = Degree_in_sample.get(Id_uvw);
							if (val >= k) {
								glob_score[Id_uvw] = glob_score[Id_uvw] + 1;
							}
						}
					}

					System.out.println("size of triangles after sampling " + Degree_in_sample.size() + "vs"
							+ list_of_triangles.size());

					// after sampling:
					Set<Integer> finalList = new HashSet<Integer>();

					for (Integer Id_uvw : Degree_in_sample.keySet()) {
						double val = (double) glob_score[Id_uvw] / (double) num_samples;
						if (val > eta) {
							// System.out.println("id " + i);
							finalList.add(Id_uvw);
						}
					}

					System.out.println(finalList.size() == list_of_triangles.size());

					BitSet visited_final = new BitSet(num_tria_alive);
					for (Integer id_uvw : finalList) {
						if (visited_final.get(id_uvw) == false) {

							List<Integer> tri_list = new ArrayList<Integer>();

							Stack<Integer> stack = new Stack<Integer>();
							stack.push(id_uvw);

							while (!(stack.isEmpty())) {

								int tri_id = stack.peek();
								stack.pop();

								int u = first_vert[tri_id];
								int v = second_vert[tri_id];
								int w = third_vert[tri_id];

								if (visited_final.get(tri_id) == false) {

									visited_final.set(tri_id);
									tri_list.add(tri_id);

									int[] u_neighbors = G.successorArray(u);
									int[] v_neighbors = G.successorArray(v);
									int[] w_neighbors = G.successorArray(w);
									int u_deg = G.outdegree(u);
									int v_deg = G.outdegree(v);
									int w_deg = G.outdegree(w);

									List<Integer> _4cli_list = Find_4clique(u, v, w, u_neighbors, v_neighbors,
											w_neighbors, u_deg, v_deg, w_deg);
									int len_list = _4cli_list.size();
									List<Integer> neigh = new ArrayList<Integer>(len_list * 3);

									for (int j = 0; j < len_list; j++) {
										int z = _4cli_list.get(j);

										int id_uvz = complex_triangle_Id(u, v, z);
										int id_uwz = complex_triangle_Id(u, w, z);
										int id_vwz = complex_triangle_Id(v, w, z);

										if ((id_uvz != -1) && (id_uwz != -1) && (id_vwz != -1)) {
											if ((finalList.contains(id_uvz)) && (finalList.contains(id_uwz))
													&& (finalList.contains(id_vwz))) {
												neigh.add(id_uvz);
												neigh.add(id_uwz);
												neigh.add(id_vwz);
											}
										}

									} // we find all the neighbors

									for (int i = 0; i < neigh.size(); i++) {
										int id_neigh = neigh.get(i);
										if (visited_final.get(id_neigh) == false) {

											stack.push(id_neigh);
										}
									}
								}
							}
							if (tri_list.size() > 1) {
								Set<Integer> edge_final = new HashSet<Integer>();
								Set<Integer> vert_final = new HashSet<Integer>();

								for (int j = 0; j < tri_list.size(); j++) {
									int id_can = tri_list.get(j);

									int u = first_vert[id_can];
									int v = second_vert[id_can];
									int w = third_vert[id_can];

									int id_uv = Arrays.binarySearch(edgeTail, offset[Math.min(u, v)],
											offset[Math.min(u, v) + 1], Math.max(u, v));
									int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u, w)],
											offset[Math.min(u, w) + 1], Math.max(u, w));
									int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(v, w)],
											offset[Math.min(v, w) + 1], Math.max(v, w));

									vert_final.add(u);
									vert_final.add(v);
									vert_final.add(w);

									edge_final.add(id_uv);
									edge_final.add(id_uw);
									edge_final.add(id_vw);
								}

								//double density = find_density(edge_final, vert_final.size());
								//double pcc = find_clust_coef(vert_final,edge_final, tri_list);
								//wg.write(density + "\t" + pcc + "\t" +  edge_final.size() + "\t" + vert_final.size() + "\n");
								for (Integer e : edge_final) {
									int u = findHead(e);
									int v = edgeTail[e];
									wg.write("\""+ u + "\"" + "->" + "\""+ v + "\"");
									wg.write("\n");
								}
								
								 //wg.write("\n");
								 wg.write("****************************************************************" +"\n");
							}

						}

						// print them as connected components
					}
				}
			} // end of the process for connected components

			System.out.println("number of connected componenets for k=" + k + " is " + count_con);
		}
		wg.close();
		//l.close();
	}

	public double find_density(Set<Integer> edge_list, int num_vertices) {
		double sum = 0.0;
		for (Integer id_e : edge_list) {
			sum += edge_prob[id_e.intValue()];
		}

		// int num_vertices = vert_list.size();
		return sum / ((double) 0.5 * num_vertices * (num_vertices - 1));
	}

	public List<Integer> intersect_3_on_samllGraph(int u, int v, int w, ArrayList<Integer> u_neighbors,
			ArrayList<Integer> v_neighbors, ArrayList<Integer> w_neighbors, int u_deg, int v_deg, int w_deg,
			int id_uvw) {

		int max_cliques_contained = Math.max(u_deg, Math.max(v_deg, w_deg));
		List<Integer> clique_list = new ArrayList<>(max_cliques_contained);

		for (int i = 0, j = 0, k = 0; i < u_deg && j < v_deg && k < w_deg;) {

			int u_neigh = u_neighbors.get(i);
			int v_neigh = v_neighbors.get(j);
			int w_neigh = w_neighbors.get(k);

			if (u_neigh == v_neigh && u_neigh == w_neigh) { // Find a 4-clique !
				int z = u_neighbors.get(i); // a 4-clique is found

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

					if ((edge_in_sample.get(id_uz) == true) && (edge_in_sample.get(id_vz) == true)
							&& (edge_in_sample.get(id_wz) == true)) {
						clique_list.add(z);
					}
				}

				i++;
				j++;
				k++;
				continue;
			}

			if (u_neigh < w_neigh) {
				if (u_neigh == v_neigh) {
					i++;
					j++;
					continue;
				} else if (u_neigh < v_neigh) {
					i++;
					continue;
				} else { // u_neighbors[i] > v_neighbors[j]
					j++;
					continue;
				}
			}

			if (u_neigh > w_neigh) {
				if (w_neigh == v_neigh) {
					k++;
					j++;
					continue;
				} else if (w_neigh < v_neigh) {
					k++;
					continue;
				} else { // w_neighbors[k] > v_neighbors[j]
					j++;
					continue;
				}
			}

			if (u_neigh == w_neigh) {
				if (u_neigh < v_neigh) {
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

	List<Integer> intersection_on_smallGraph(int u, int v, List<Integer> u_neighbors, List<Integer> v_neighbors,
			int u_deg, int v_deg) {

		List<Integer> wlist = new ArrayList<Integer>();
		for (int i = 0, j = 0; i < u_deg && j < v_deg;) {

			int u_com = u_neighbors.get(i);
			int v_com = v_neighbors.get(j);

			if (u_com == v_com) { // Find a triangle !
				// the if condition below is to avoid double counting:
				if (u_com > v) {

					int w = u_neighbors.get(i);
					wlist.add(w);

				}
				i++;
				j++;
				continue;
			}

			if (u_com < v_com) {
				i++;
				continue;
			}

			if (u_com > v_com) {
				j++;
				continue;
			}
		}

		return wlist;
	}

	/*
	 * public void final_result() {
	 * 
	 * double sum_dens = 0.0; double max_dens = 0.0; for (int i = 0; i <
	 * density.size(); i++) { double dens_value = density.get(i);
	 * 
	 * if (dens_value > max_dens) { max_dens = dens_value; }
	 * 
	 * sum_dens += dens_value; }
	 * 
	 * double ave_dens = sum_dens / ((double) density.size());
	 * 
	 * System.out.println("average density is: " + ave_dens);
	 * System.out.println("maximum density is: " + max_dens);
	 * 
	 * double sum_cc = 0.0; double max_cc = 0.0; for (int i = 0; i <
	 * clust_coeif.size(); i++) { double cc_value = clust_coeif.get(i);
	 * 
	 * if (cc_value > max_cc) { max_cc = cc_value; }
	 * 
	 * sum_cc += cc_value; }
	 * 
	 * double ave_cc = sum_cc / ((double) clust_coeif.size());
	 * 
	 * System.out.println("average clustering coefficient is: " + ave_cc);
	 * System.out.println("maximum clustering coefficient is: " + max_cc);
	 * 
	 * double sum_vertices = 0.0; int max_vertices = 0;
	 * 
	 * for (int i = 0; i < num_of_vertices.size(); i++) { int vertices_value =
	 * num_of_vertices.get(i);
	 * 
	 * if (vertices_value > max_vertices) { max_vertices = vertices_value; }
	 * 
	 * sum_vertices += vertices_value; }
	 * 
	 * double ave_vertices = sum_vertices / ((double) num_of_vertices.size());
	 * 
	 * System.out.println("average number of vertices is: " + ave_vertices);
	 * System.out.println("maximum number of vertices is: " + max_vertices);
	 * 
	 * double sum_edges = 0.0; int max_edges = 0;
	 * 
	 * for (int i = 0; i < num_of_edges.size(); i++) { int edges_value =
	 * num_of_edges.get(i);
	 * 
	 * if (edges_value > max_edges) { max_edges = edges_value; }
	 * 
	 * sum_edges += edges_value; }
	 * 
	 * double ave_edges = sum_edges / ((double) num_of_edges.size());
	 * 
	 * System.out.println("average number of edges is: " + ave_edges);
	 * System.out.println("maximum number of edges is: " + max_edges); }
	 */

	public double find_density(BitSet edge_list, Set<Integer> vert_list) {
		double sum = 0.0;
		for (int e = 0; e < m; e++) {
			if (edge_list.get(e) == true) {
				sum += edge_prob[e];
			}
		}
		int num_vertices = vert_list.size();
		return sum / ((double) 0.5 * num_vertices * (num_vertices - 1));
	}

	public void DFS_iter(int id, boolean[] visited, Set<Integer> s, List<Integer> results) {

		Stack<Integer> stack = new Stack<Integer>();
		stack.push(id);

		while (stack.empty() == false) {

			id = stack.peek();
			stack.pop();

			int u = first_vert[id];
			int v = second_vert[id];
			int w = third_vert[id];

			if (visited[id] == false) {

				visited[id] = true;
				// System.out.print(" " + id + " ");
				results.add(id);
			}

			int[] u_neighbors = G.successorArray(u);
			int[] v_neighbors = G.successorArray(v);
			int[] w_neighbors = G.successorArray(w);
			int u_deg = G.outdegree(u);
			int v_deg = G.outdegree(v);
			int w_deg = G.outdegree(w);

			List<Integer> _4cli_list = Find_4clique(u, v, w, u_neighbors, v_neighbors, w_neighbors, u_deg, v_deg,
					w_deg);
			int len = _4cli_list.size();
			List<Integer> neigh = new ArrayList<Integer>(len * 3);

			for (int j = 0; j < len; j++) {
				int z = _4cli_list.get(j);

				int id_uvz = complex_triangle_Id(u, v, z);
				int id_uwz = complex_triangle_Id(u, w, z);
				int id_vwz = complex_triangle_Id(v, w, z);

				if ((id_uvz != -1) && (id_uwz != -1) && (id_vwz != -1)) {
					if ((s.contains(id_uvz)) && (s.contains(id_uwz)) && (s.contains(id_vwz))) {
						neigh.add(id_uvz);
						neigh.add(id_uwz);
						neigh.add(id_vwz);
					}
				}

			} // we find all the neighbors

			for (int i = 0; i < neigh.size(); i++) {
				int id_neigh = neigh.get(i);
				if (visited[id_neigh] == false) {

					stack.push(id_neigh);
				}
			}

		}
	}

	public void sampling() {
		for (int e = 0; e < m; e++) {
			if (member_of_component_for_edge.get(e) == true) { // if the edge is part of k-nucleus graph
				double prob_e = edge_prob[e];
				long random_num = (long) (Math.pow(10, precision) * Math.random());
				double random_prob = (double) random_num * Math.pow(10, -precision);
				// System.out.println(random_prob);
				if ((random_prob) < prob_e) { // we pick the edge
					edge_in_sample.set(e);
				}
			}
		}
	}

	List<Integer> Find_4clique(int u, int v, int w, int[] u_neighbors, int[] v_neighbors, int[] w_neighbors, int u_deg,
			int v_deg, int w_deg) {

		int minlength = 0;
		minlength = Math.min(u_deg, Math.min(v_deg, w_deg));
		List<Integer> clique_list = new ArrayList<Integer>(minlength);

		for (int i = 0, j = 0, k = 0; i < u_deg && j < v_deg && k < w_deg;) {

			if (u_neighbors[i] == v_neighbors[j] && u_neighbors[i] == w_neighbors[k]) { // Find a 4-clique !
//            System.out.println(u + ", " + v + ", " + w  + ", " + u_neighbors[i]); // for checking
				clique_list.add(u_neighbors[i]);
				// degree++;
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

	public int ReadFile(String filename, int k_value) throws Exception {

		BufferedReader br = new BufferedReader(new FileReader(filename));
		String line = br.readLine();
		int index_id = 0;
		int num = 0;

		while (line != null) {

			int kValue = Integer.parseInt(line.split("\t")[3]);
			if (kValue >= k_value) {

				num++;
				// candidate.add(index_id);

				int u = Integer.parseInt(line.split("\t")[0]);
				int v = Integer.parseInt(line.split("\t")[1]);
				int w = Integer.parseInt(line.split("\t")[2]);

				int id_uv = Arrays.binarySearch(edgeTail, offset[Math.min(u, v)], offset[Math.min(u, v) + 1],
						Math.max(u, v));
				int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u, w)], offset[Math.min(u, w) + 1],
						Math.max(u, w));
				int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(v, w)], offset[Math.min(v, w) + 1],
						Math.max(v, w));

				// edge_in_knucleus.set(id_uv);
				// edge_in_knucleus.set(id_uw);
				// edge_in_knucleus.set(id_vw);

			}

			index_id++;
			line = br.readLine();
		}

		br.close();
		return num;
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

	public void Populate_edges() {
		int startPos = 0;
		for (int u = 0; u < n; u++) {
			offset[u] = startPos;

			int[] u_successor = G.successorArray(u);
			Label[] u_label = G.labelArray(u);

			int count = 0;
			for (int j = 0; j < u_successor.length; j++) {
				int v = u_successor[j];
				if (u < v) {
					edgeTail[startPos + count] = v;
					edge_prob[startPos + count] = u_label[j].getLong() * (1.0 / Math.pow(10, precision));
					count++;
				}
			}
			startPos += count;
		}
		offset[n] = offset[n - 1];
	}

	List<Integer> intersection_on_SmallGraph(int u, int v, int id_uv, List<Integer> neighborlist_u,
			List<Integer> neighborlist_v, int u_deg, int v_deg) {

		List<Integer> wlist = new ArrayList<Integer>();
		for (int i = 0, j = 0; i < u_deg && j < v_deg;) {

			int u_neigh = neighborlist_u.get(i);
			int v_neigh = neighborlist_v.get(j);

			if (u_neigh == v_neigh) { // Find a triangle !
				// the if condition below is to avoid double counting:
				if (u_neigh > v) {

					int w = neighborlist_u.get(i);

					int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u, w)], offset[Math.min(u, w) + 1],
							Math.max(u, w));
					int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(w, v)], offset[Math.min(w, v) + 1],
							Math.max(w, v));
					double exiprob = edge_prob[id_uv] * edge_prob[id_uw] * edge_prob[id_vw];

					if (exiprob > eta) {

						if ((edge_in_sample.get(id_uw) == true) && (edge_in_sample.get(id_vw) == true)) {
							wlist.add(w);
						}
					}
				}
				i++;
				j++;
				continue;
			}

			if (u_neigh < v_neigh) {
				i++;
				continue;
			}

			if (u_neigh > v_neigh) {
				j++;
				continue;
			}
		}

		return wlist;
	}

	List<Integer> intersection(int u, int v, int id_uv, int[] u_neighbors, int[] v_neighbors, int u_deg, int v_deg) {

		int index = 0;
		List<Integer> wlist = new ArrayList<Integer>();
		for (int i = 0, j = 0; i < u_deg && j < v_deg;) {

			if (u_neighbors[i] == v_neighbors[j]) { // Find a triangle !
				// the if condition below is to avoid double counting:
				if (u_neighbors[i] > v) {

					int w = u_neighbors[i];

					int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u, w)], offset[Math.min(u, w) + 1],
							Math.max(u, w));
					int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(w, v)], offset[Math.min(w, v) + 1],
							Math.max(w, v));
					double exiprob = edge_prob[id_uv] * edge_prob[id_uw] * edge_prob[id_vw];

					if (exiprob > eta) {

						if ((edge_in_sample.get(id_uw) == true) && (edge_in_sample.get(id_vw) == true)) {
							wlist.add(u_neighbors[i]);
							index++;
						}
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

	public static void main(String[] args) throws Exception {
		
		String basename = args[0];
		//long tircount = Long.valueOf(args[1]);
		int precision = Integer.valueOf(args[1]);
		double eta = Double.valueOf(args[2]);
		int num_sample = Integer.valueOf(args[3]);

		//String basename = "krogan-proc.w";
		//int num_sample = 200;


		//String basename = "krogan-proc.w";

		//double eta = 0.1;
		//int num_sample = 2000;
		// This is the final version

		// String basename = "DBLP_algorithm-proc.w";
		//int precision = 3;
		
		WeaklyGlobalNucleus GNC = new WeaklyGlobalNucleus(basename,
				precision, eta, num_sample);
		GNC.Main_BasicSampling();

	}

}
