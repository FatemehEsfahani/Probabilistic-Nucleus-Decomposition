import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.SortedMap;
import java.util.Stack;
import java.util.TreeMap;
import java.util.stream.IntStream;

import it.unimi.dsi.webgraph.ImmutableGraph;
import it.unimi.dsi.webgraph.labelling.ArcLabelledImmutableGraph;
import it.unimi.dsi.webgraph.labelling.Label;

//This code I have written for task-driven team formation.
//1) It reads the nucleus scores obtained by DP.
//2) It finds the connected components.
//3) for each connected component which contains the authors of interest, it finds the global nuclei. In this version, I store the 
//connected component as a new graph.

public class Golbal_nuclei_Finding_onSmallGraph_ForExperiments {

	String basename;
	int precision;
	double eta;
	
	//double eta2;

	ArcLabelledImmutableGraph G; // uncertain graph
	int n; // vertices of G
	int m; // edges of G

	String filename; // nucleus scores of triangles
	int num_samples;

	final int[] edgeTail;
	final int[] offset; // offset array for edges
	double[] edge_prob;
	
	//String globalfile;
	//String localfile;

	//String file = "globalResult2.txt";
	//BufferedWriter wr;

	// for each triangle (u,v,w), where u<v<w:
	int[] first_vert; // first_vert stores u
	int[] second_vert; // second_vert stores v
	int[] third_vert; // third_vert stores
	int[] offset_uv; // offset array for triangles
	Double[] triangle_prob;

	double epsilon0 = 1.0E-16;
	double epsilon1 = 1.0 - epsilon0;

	BigInteger t_total_cnt = BigInteger.ZERO;
	int num_tria_alive;
	
	int max_score;
	//String keyword;

	HashMap<Integer, List<Integer>> triangle_comneigh = new HashMap<>();
	HashMap<Integer, List<Integer>> edge_comneigh = new HashMap<>();

	// *********for connected component finding:
	int count_CC = 0; // number of connected component
	
	BitSet visited;

	int[] support; // stores 4kdegree for each triangle

	List<Integer> list_of_vertices;
	List<Integer> list_of_triangles;
	List<Integer> list_of_edges;

	BitSet member_of_component;
	BitSet member_of_component_for_edge;
	
	
	//String[] author_id;

	

	ArrayList<BitSet> edge_in_sample; // edges alive in each sample graph
	// int[] sortedTrinagle; // stores triangles in ascending order of their
	// 4kdegrees

	// ArrayList<List<Integer>> triangle_clique;
	// ArrayList<List<Integer>> edge_triangle;

	
	//public Golbal_nuclei_Finding_onSmallGraph_ForExperiments(String basename, int precision, double eta, int num_samples,double eta2) throws Exception {
	public Golbal_nuclei_Finding_onSmallGraph_ForExperiments(String basename, int precision, double eta, int num_samples) throws Exception {

		this.basename = basename;
		this.precision = precision;
		this.eta = eta;
		//this.eta2 = eta2;
		
		this.filename = this.basename + "-" + this.eta + "-" + "finalsupp_DP_Modified.txt";
		//this.keyword = keyword;
		
		G = ArcLabelledImmutableGraph.load(basename);
		n = G.numNodes();
		m = (int) (G.numArcs() / 2);

		edgeTail = new int[m]; // store tail for each edge
		offset = new int[n + 1]; // edge offset for edgeHead array
		edge_prob = new double[m];

		this.num_samples = num_samples;
		
		//author_id = new String[n];

		System.out.println("start storing edges...");
		Populate_edges(); // Store edges

		System.out.println("start enumerating triangles...");
		EnumerateTriangles(); // find triangles with prob > eta

		edge_in_sample = new ArrayList<BitSet>(this.num_samples);
		for (int N = 0; N < this.num_samples; N++) {
			BitSet bitarray = new BitSet(m);
			edge_in_sample.add(N, bitarray);
		}

		System.out.println("start sampling...");
		sampling();

		offset_uv = new int[m + 1]; // starting index of triangles with edge (u,v) where u<v
		first_vert = new int[num_tria_alive];
		second_vert = new int[num_tria_alive];
		third_vert = new int[num_tria_alive];
		triangle_prob = new Double[num_tria_alive];

		System.out.println("start storing triangles...");
		Populate_triangles();

		visited = new BitSet(num_tria_alive);
		support = new int[num_tria_alive];

		list_of_vertices = new ArrayList<>();
		list_of_triangles = new ArrayList<>();
		list_of_edges = new ArrayList<>();

		member_of_component = new BitSet(n);
		member_of_component_for_edge = new BitSet(m);

		System.out.println("start reading nucleus scores...");
		ReadFile(filename);
		
		int max=0;
		for (int i=0;i<support.length;i++) {
			int val = support[i];
			if (val>= max) {
				max = val;
			}
		}
		max_score = max;
		
		System.out.println("maximum nucleuss score is: " + max_score);

	}

	public void sampling() {
		
		double minWeight = 1;
		double maxWeight = Math.pow(10, precision);
		
		for (int N = 0; N < num_samples; N++) {

			for (int e = 0; e < m; e++) {
				
				long random_num = (long)(minWeight + (maxWeight-minWeight)*Math.random());

				double prob_e = edge_prob[e];
				
				//long random_num = (long) (Math.pow(10, precision) * Math.random());

				//long random_num = (long) (100 * Math.random());
				//double random_prob = (double) random_num * 0.01;
				
				double random_prob = (double) random_num * Math.pow(10, -precision);

				if ((random_prob) < prob_e) { // we pick the edge
					// edge_in_sample.get(e).set(N);
					edge_in_sample.get(N).set(e);
				}

			}
		}
	} // end sampling method

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

	public void findConnectedComponnents() throws Exception {
		
		String globalfile = "Global-"+num_samples + "-" +basename + "-" + eta + "-.txt"; 
		BufferedWriter wg = new BufferedWriter(new FileWriter(globalfile));
		
		System.out.println("start finding connected components and their global subgraph....");
		
		long startTime = System.currentTimeMillis();
		for (int k=max_score;k>0;k--) {
			
			//if (k==1) {
			//	break;
			//}
			
			System.out.println("k-nucleus=" + k);
			
			
			//wg.write("k-nucleus " + k + "=" + "\n");
			int count_con = 0;
			visited.clear();
			
			for (int id = 0; id < num_tria_alive; id++) {

				if ((visited.get(id) == false) && (support[id] >= k)) {

					list_of_vertices = new ArrayList<>();
					list_of_triangles = new ArrayList<>();
					list_of_edges = new ArrayList<>();
					
					member_of_component.clear();
					member_of_component_for_edge.clear();

					DFS_iter(id,k);
					// now one connected component has been detected
					count_con++;
					
					if (count_con >= 1000) {
						break; //enough connected component for global has been detected.
					}
				
					
					//if (authors_in_comp == true) {
						
						//Perform global decomposition
						// each time we make it empty for new components
						triangle_comneigh.clear();
						edge_comneigh.clear();

						// create induced graph based on vertcies in the componnent
						SortedMap<Integer, SortedMap<Integer, Double>> G_conn = new TreeMap<>();
						for (int i = 0; i < list_of_edges.size(); i++) {
							int id_e = list_of_edges.get(i);
							int u = findHead(id_e);
							int v = edgeTail[id_e];
							
							//store connected component
							//wl.write("(" + author_id[u] + "," + author_id[v] + ")" + " ");
							//wl.write("\""+ author_id[u] + "\"" + "->" + "\""+ author_id[v] + "\"");
							//wl.write("\n");

							if (!G_conn.containsKey(u)) {
								G_conn.put(u, new TreeMap<>());
							}

							if (!G_conn.containsKey(v)) {
								G_conn.put(v, new TreeMap<>());
							}

							G_conn.get(u).put(v, edge_prob[id_e]);
							G_conn.get(v).put(u, edge_prob[id_e]);

						}
						
						//wl.write("\n");
						//wl.write("****************************************************************" + "\n");

						for (Integer u : G_conn.keySet()) {
							for (Integer v : G_conn.get(u).keySet()) {
								if (u.intValue() < v.intValue()) {
									ArrayList<Integer> u_neigh = new ArrayList<Integer>(G_conn.get(u).keySet());
									ArrayList<Integer> v_neigh = new ArrayList<Integer>(G_conn.get(v).keySet());

									int u_deg = u_neigh.size();
									int v_deg = v_neigh.size();
									List<Integer> tri_lits = intersection_on_smallGraph(u.intValue(), v.intValue(), u_neigh,
											v_neigh, u_deg, v_deg);

									int index = tri_lits.size();
									for (int j = 0; j < index; j++) {
										int w = tri_lits.get(j);
										int id_uvw = findTrinagle_Id(u, v, w);
										ArrayList<Integer> w_neigh = new ArrayList<Integer>(G_conn.get(w).keySet());
										int w_deg = w_neigh.size();
										List<Integer> arr = intersect_3_on_samllGraph(u.intValue(), v.intValue(), w, u_neigh,
												v_neigh, w_neigh, u_deg, v_deg, w_deg, id_uvw);
										triangle_comneigh.put(id_uvw, arr);

									}

								}
							}
						}

						for (int i = 0; i < list_of_edges.size(); i++) {
							int id_e = list_of_edges.get(i);
							int u = findHead(id_e);
							int v = edgeTail[id_e];
							ArrayList<Integer> u_neigh = new ArrayList<Integer>(G_conn.get(u).keySet());
							ArrayList<Integer> v_neigh = new ArrayList<Integer>(G_conn.get(v).keySet());

							int u_deg = u_neigh.size();
							int v_deg = v_neigh.size();

							List<Integer> arr = FindCommNeigh_for_edge_onSmallGraph(u, v, u_deg, v_deg, u_neigh, v_neigh);
							edge_comneigh.put(id_e, arr);
						}
						
						Main_globalDecomposition(wg,k);

					//}
					

				}

			}
			long endTime = System.currentTimeMillis();
			System.out.println("Time elapsed (sec) for loading the graph = " + (endTime - startTime) / 1000.0);
			System.out.println("number of connected componenets for k=" + k + " is " + count_con);
		}
		
		wg.close();
		//wl.close();
		

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

	public List<Integer> FindCommNeigh_for_edge_onSmallGraph(int u, int v, int u_deg, int v_deg, List<Integer> u_neigh,
			List<Integer> v_neigh) {

		List<Integer> li = new ArrayList<Integer>();

		for (int i = 0, j = 0; i < u_deg && j < v_deg;) {
			int u_N = u_neigh.get(i);
			int v_N = v_neigh.get(j);

			if (u_N == v_N) { // Find a triangle !

				int w = u_neigh.get(i); // triangle (u,v,w)

				int id_uvw = complex_triangle_Id(u, v, w);
				if (id_uvw != -1) {

					li.add(id_uvw);
				}
				i++;
				j++;
				continue;
			}

			if (u_N < v_N) {
				i++;
				continue;
			}

			if (u_N > v_N) {
				j++;
				continue;
			}
		}

		return li;
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

					if ((member_of_component_for_edge.get(id_uz) == true)
							&& (member_of_component_for_edge.get(id_vz) == true)
							&& (member_of_component_for_edge.get(id_wz) == true)) {
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

	public void Main_globalDecomposition(BufferedWriter wg, int k_value) throws IOException {

		
		BitSet member_in_C_k = new BitSet(num_tria_alive);
		for (int i = 0; i < list_of_triangles.size(); i++) {
			int id_uvw = list_of_triangles.get(i);
			member_in_C_k.set(id_uvw);
		}

		// edges alive in C_k
		// BitSet edge_in_C_k = new BitSet(m);
		// member_of_component_for_edge instead of edge_in_C_k

		//System.out.println("number of edges are: " + member_of_component_for_edge.cardinality());
		//System.out.println("number of candidate triangle are: " + member_in_C_k.cardinality());
		//System.out.println(
			//	"this value should be equal to the length C_k array list which is :" + list_of_triangles.size());

		System.out.println("start checking the degree contraint for each candidate traingle...");
		long startTime_DC = System.currentTimeMillis();
		degree_constraint_check(list_of_triangles, member_in_C_k, member_of_component_for_edge,k_value);
		long endTime_DC = System.currentTimeMillis();
		System.out.println("Time elapsed (sec) for creating candidates = " + (endTime_DC - startTime_DC) / 1000.0);
		// Now, each triangle has enough degree (support)

		System.out.println("start creating the graph Q....");
		long startTime_CGR = System.currentTimeMillis();
		if (member_in_C_k.cardinality() != 0) { // if there exists any candidate

			System.out.println("processing " + list_of_triangles.size() + " of triangles.");

			for (int j = 0; j < list_of_triangles.size(); j++) {
				if (j % 100 == 0) {
					System.out.println(j);
				}
				int id = list_of_triangles.get(j);
				if (member_in_C_k.get(id) == true) {

					boolean hold_condition = true;

					Queue<Integer> q = new LinkedList<Integer>();
					BitSet edge_in_Q = new BitSet(m);
					Set<Integer> triangle_in_Q = new HashSet<Integer>(num_tria_alive);
					
					boolean[] visited = new boolean[num_tria_alive];
					triangle_in_Q.add(id);
					q.add(id);
					visited[id] = true;
					while (!(q.isEmpty())) {
						int id_uvw = q.remove();
						// triangle_in_Q.add(id_uvw);

						int u = first_vert[id_uvw];
						int v = second_vert[id_uvw];
						int w = third_vert[id_uvw];

						int id_uv = Arrays.binarySearch(edgeTail, offset[Math.min(u, v)], offset[Math.min(u, v) + 1],
								Math.max(u, v));

						if (edge_in_Q.get(id_uv) == false) {
							edge_in_Q.set(id_uv);

							List<Integer> tri_list = edge_comneigh.get(id_uv);
							for (int ii = 0; ii < tri_list.size(); ii++) {

								// int com_neigh = tri_list.get(ii);
								int triangle_id = tri_list.get(ii);

								if (triangle_id != id_uvw) {
									if (visited[triangle_id] == false) {
										int com_neigh;
										if ((first_vert[triangle_id] != u) && (first_vert[triangle_id] != v)) {
											com_neigh = first_vert[triangle_id];
										} else {
											if (second_vert[triangle_id] == v) {
												com_neigh = third_vert[triangle_id];
											} else {
												com_neigh = second_vert[triangle_id];
											}
										}

										int id_ucom = Arrays.binarySearch(edgeTail, offset[Math.min(u, com_neigh)],
												offset[Math.min(u, com_neigh) + 1], Math.max(u, com_neigh));
										int id_vcom = Arrays.binarySearch(edgeTail, offset[Math.min(v, com_neigh)],
												offset[Math.min(v, com_neigh) + 1], Math.max(v, com_neigh));

										if ((edge_in_Q.get(id_ucom) == true) && (edge_in_Q.get(id_vcom) == true)) {
											triangle_in_Q.add(triangle_id);
											q.add(triangle_id);
											visited[triangle_id] = true;
										}
									}
								}
							}
						}

						int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u, w)], offset[Math.min(u, w) + 1],
								Math.max(u, w));

						if (edge_in_Q.get(id_uw) == false) {
							edge_in_Q.set(id_uw);
							List<Integer> tri_list = edge_comneigh.get(id_uw);
							for (int ii = 0; ii < tri_list.size(); ii++) {
								// int com_neigh = tri_list.get(ii);

								int triangle_id = tri_list.get(ii);
								if (triangle_id != id_uvw) {
									if (visited[triangle_id] == false) {
										int com_neigh;
										if ((first_vert[triangle_id] != u) && (first_vert[triangle_id] != w)) {
											com_neigh = first_vert[triangle_id];
										} else {
											if (second_vert[triangle_id] == w) {
												com_neigh = third_vert[triangle_id];
											} else {
												com_neigh = second_vert[triangle_id];
											}
										}

										int id_ucom = Arrays.binarySearch(edgeTail, offset[Math.min(u, com_neigh)],
												offset[Math.min(u, com_neigh) + 1], Math.max(u, com_neigh));
										int id_wcom = Arrays.binarySearch(edgeTail, offset[Math.min(w, com_neigh)],
												offset[Math.min(w, com_neigh) + 1], Math.max(w, com_neigh));
										if ((edge_in_Q.get(id_ucom) == true) && (edge_in_Q.get(id_wcom) == true)) {
											triangle_in_Q.add(triangle_id);
											q.add(triangle_id);
											visited[triangle_id] = true;
										}

									}
								}

							}
						}

						int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(v, w)], offset[Math.min(v, w) + 1],
								Math.max(v, w));
						if (edge_in_Q.get(id_vw) == false) {
							edge_in_Q.set(id_vw);
							List<Integer> tri_list = edge_comneigh.get(id_vw);
							for (int ii = 0; ii < tri_list.size(); ii++) {
								// int com_neigh = tri_list.get(ii);
								int triangle_id = tri_list.get(ii);

								if (triangle_id != id_uvw) {
									if (visited[triangle_id] == false) {
										int com_neigh;

										if ((first_vert[triangle_id] != v) && (first_vert[triangle_id] != w)) {
											com_neigh = first_vert[triangle_id];
										} else {
											if (second_vert[triangle_id] == w) {
												com_neigh = third_vert[triangle_id];
											} else {
												com_neigh = second_vert[triangle_id];
											}
										}

										int id_vcom = Arrays.binarySearch(edgeTail, offset[Math.min(v, com_neigh)],
												offset[Math.min(v, com_neigh) + 1], Math.max(v, com_neigh));
										int id_wcom = Arrays.binarySearch(edgeTail, offset[Math.min(w, com_neigh)],
												offset[Math.min(w, com_neigh) + 1], Math.max(w, com_neigh));
										if ((edge_in_Q.get(id_vcom) == true) && (edge_in_Q.get(id_wcom) == true)) {
											triangle_in_Q.add(triangle_id);
											q.add(triangle_id);
											visited[triangle_id] = true;
										}

									}

								}

							}
						}

						// third version
						List<Integer> cli = triangle_comneigh.get(id_uvw);
						
						int len = cli.size();
						int count = 0;
						int jj = 0;

						while ((jj < len) && (count < k_value)) {
							int z = cli.get(jj);
							int id_uz = Arrays.binarySearch(edgeTail, offset[Math.min(u, z)],
									offset[Math.min(u, z) + 1], Math.max(u, z));
							int id_vz = Arrays.binarySearch(edgeTail, offset[Math.min(v, z)],
									offset[Math.min(v, z) + 1], Math.max(v, z));
							int id_wz = Arrays.binarySearch(edgeTail, offset[Math.min(w, z)],
									offset[Math.min(w, z) + 1], Math.max(w, z));
							if ((edge_in_Q.get(id_uz) == true) && (edge_in_Q.get(id_vz) == true)
									&& (edge_in_Q.get(id_wz) == true)) {
								// the clique already exists
								count++;
							}
							jj++;
						}
						jj = 0;

						while ((jj < len) && (count < k_value)) {
							int z = cli.get(jj);
							int id_uz = Arrays.binarySearch(edgeTail, offset[Math.min(u, z)],
									offset[Math.min(u, z) + 1], Math.max(u, z));
							int id_vz = Arrays.binarySearch(edgeTail, offset[Math.min(v, z)],
									offset[Math.min(v, z) + 1], Math.max(v, z));
							int id_wz = Arrays.binarySearch(edgeTail, offset[Math.min(w, z)],
									offset[Math.min(w, z) + 1], Math.max(w, z));
							if ((member_of_component_for_edge.get(id_uz) == true)
									&& (member_of_component_for_edge.get(id_vz) == true)
									&& (member_of_component_for_edge.get(id_wz) == true)) {
								if ((edge_in_Q.get(id_uz) == true) && (edge_in_Q.get(id_vz) == true)
										&& (edge_in_Q.get(id_wz) == true)) {

								} else {
									count++;
								}
								int id_uvz = complex_triangle_Id(u, v, z);
								int id_uwz = complex_triangle_Id(u, w, z);
								int id_vwz = complex_triangle_Id(v, w, z);

								if (visited[id_uvz] == false) {
									triangle_in_Q.add(id_uvz);
									q.add(id_uvz);
									visited[id_uvz] = true;
								}
								if (visited[id_uwz] == false) {
									triangle_in_Q.add(id_uwz);
									q.add(id_uwz);
									visited[id_uwz] = true;
								}
								if (visited[id_vwz] == false) {
									triangle_in_Q.add(id_vwz);
									q.add(id_vwz);
									visited[id_vwz] = true;
								}

								if (edge_in_Q.get(id_uz) == false) {
									edge_in_Q.set(id_uz);
									List<Integer> li_tri = edge_comneigh.get(id_uz);
									for (int dd = 0; dd < li_tri.size(); dd++) {

										int triangle_id = li_tri.get(dd);
										if ((triangle_id != id_uvz) && (triangle_id != id_uwz)) {
											if (visited[triangle_id] == false) {
												int u_1 = first_vert[triangle_id];
												int v_1 = second_vert[triangle_id];
												int w_1 = third_vert[triangle_id];

												int id_u1v1 = Arrays.binarySearch(edgeTail, offset[Math.min(u_1, v_1)],
														offset[Math.min(u_1, v_1) + 1], Math.max(u_1, v_1));
												int id_u1w1 = Arrays.binarySearch(edgeTail, offset[Math.min(u_1, w_1)],
														offset[Math.min(u_1, w_1) + 1], Math.max(u_1, w_1));
												int id_v1w1 = Arrays.binarySearch(edgeTail, offset[Math.min(v_1, w_1)],
														offset[Math.min(v_1, w_1) + 1], Math.max(v_1, w_1));

												if ((id_uz != id_u1v1) & (id_uz != id_u1w1)) { // id_uz = id_v1w1
													if ((edge_in_Q.get(id_u1v1) == true)
															&& (edge_in_Q.get(id_u1w1) == true)) {
														triangle_in_Q.add(triangle_id);
														q.add(triangle_id);
														visited[triangle_id] = true;
													}
												} else if ((id_uz != id_u1v1) && (id_uz != id_v1w1)) { // id_uz =
																										// id_u1w1
													if ((edge_in_Q.get(id_u1v1) == true) && (edge_in_Q.get(id_v1w1))) {
														triangle_in_Q.add(triangle_id);
														q.add(triangle_id);
														visited[triangle_id] = true;
													}
												} else { // id_uz = u1v1
													if ((edge_in_Q.get(id_u1w1) == true) && (edge_in_Q.get(id_v1w1))) {
														triangle_in_Q.add(triangle_id);
														q.add(triangle_id);
														visited[triangle_id] = true;
													}

												}
											}
										}

									}
								}

								if (edge_in_Q.get(id_vz) == false) {
									edge_in_Q.set(id_vz);
									List<Integer> li_tri = edge_comneigh.get(id_vz);
									for (int dd = 0; dd < li_tri.size(); dd++) {

										int triangle_id = li_tri.get(dd);
										if ((triangle_id != id_uvz) && (triangle_id != id_vwz)) {
											if (visited[triangle_id] == false) {
												int u_1 = first_vert[triangle_id];
												int v_1 = second_vert[triangle_id];
												int w_1 = third_vert[triangle_id];

												int id_u1v1 = Arrays.binarySearch(edgeTail, offset[Math.min(u_1, v_1)],
														offset[Math.min(u_1, v_1) + 1], Math.max(u_1, v_1));
												int id_u1w1 = Arrays.binarySearch(edgeTail, offset[Math.min(u_1, w_1)],
														offset[Math.min(u_1, w_1) + 1], Math.max(u_1, w_1));
												int id_v1w1 = Arrays.binarySearch(edgeTail, offset[Math.min(v_1, w_1)],
														offset[Math.min(v_1, w_1) + 1], Math.max(v_1, w_1));

												if ((id_vz != id_u1v1) & (id_vz != id_u1w1)) { // id_vz = id_v1w1
													if ((edge_in_Q.get(id_u1v1) == true)
															&& (edge_in_Q.get(id_u1w1) == true)) {
														triangle_in_Q.add(triangle_id);
														q.add(triangle_id);
														visited[triangle_id] = true;
													}
												} else if ((id_vz != id_u1v1) && (id_vz != id_v1w1)) { // id_vz =
																										// id_u1w1
													if ((edge_in_Q.get(id_u1v1) == true) && (edge_in_Q.get(id_v1w1))) {
														triangle_in_Q.add(triangle_id);
														q.add(triangle_id);
														visited[triangle_id] = true;
													}
												} else { // id_vz = u1v1
													if ((edge_in_Q.get(id_u1w1) == true) && (edge_in_Q.get(id_v1w1))) {
														triangle_in_Q.add(triangle_id);
														q.add(triangle_id);
														visited[triangle_id] = true;
													}

												}
											}
										}

									}
								}

								if (edge_in_Q.get(id_wz) == false) {
									edge_in_Q.set(id_wz);
									List<Integer> li_tri = edge_comneigh.get(id_wz);
									for (int dd = 0; dd < li_tri.size(); dd++) {

										int triangle_id = li_tri.get(dd);
										if ((triangle_id != id_uwz) && (triangle_id != id_vwz)) {
											if (visited[triangle_id] == false) {
												int u_1 = first_vert[triangle_id];
												int v_1 = second_vert[triangle_id];
												int w_1 = third_vert[triangle_id];

												int id_u1v1 = Arrays.binarySearch(edgeTail, offset[Math.min(u_1, v_1)],
														offset[Math.min(u_1, v_1) + 1], Math.max(u_1, v_1));
												int id_u1w1 = Arrays.binarySearch(edgeTail, offset[Math.min(u_1, w_1)],
														offset[Math.min(u_1, w_1) + 1], Math.max(u_1, w_1));
												int id_v1w1 = Arrays.binarySearch(edgeTail, offset[Math.min(v_1, w_1)],
														offset[Math.min(v_1, w_1) + 1], Math.max(v_1, w_1));

												if ((id_wz != id_u1v1) & (id_wz != id_u1w1)) { // id_wz = id_v1w1
													if ((edge_in_Q.get(id_u1v1) == true)
															&& (edge_in_Q.get(id_u1w1) == true)) {
														triangle_in_Q.add(triangle_id);
														q.add(triangle_id);
														visited[triangle_id] = true;
													}
												} else if ((id_wz != id_u1v1) && (id_wz != id_v1w1)) { // id_wz =
																										// id_u1w1
													if ((edge_in_Q.get(id_u1v1) == true) && (edge_in_Q.get(id_v1w1))) {
														triangle_in_Q.add(triangle_id);
														q.add(triangle_id);
														visited[triangle_id] = true;
													}
												} else { // id_wz = u1v1
													if ((edge_in_Q.get(id_u1w1) == true) && (edge_in_Q.get(id_v1w1))) {
														triangle_in_Q.add(triangle_id);
														q.add(triangle_id);
														visited[triangle_id] = true;
													}

												}
											}
										}

									}
								}

							}
							jj++;
						}

						if (count < k_value) {
							hold_condition = false;
							break;
						}
					}
					
					// first, check if Q is not empty.
					if ((!(triangle_in_Q.isEmpty())) && (hold_condition == true)) {
						
						Check_condition(wg, triangle_in_Q, edge_in_Q,k_value);
						

					}
				}

			} // endfor

		}

		long endTime_CGR = System.currentTimeMillis();
		System.out.println("Time elapsed (sec) for main decomposition = " + (endTime_CGR - startTime_CGR) / 1000.0);
		
	}

	public void degree_constraint_check(List<Integer> C_k, BitSet member_in_C_k, BitSet edge_in_C_k, int k_value) {

		boolean update = true;
		while (update) {
			if (member_in_C_k.cardinality() == 0) {
				break;
			} else {
				update = false; // later change it if you see it should be executed earlier!
				for (int i = 0; i < C_k.size(); i++) {
					int id_uvw = C_k.get(i);
					if (member_in_C_k.get(id_uvw) == true) {

						// first, we make sure that the removal of other triangles does not cause
						// the removal of this triangle.
						// To do so, we should check whether edges in those triangles are still alive.
						int u = first_vert[id_uvw];
						int v = second_vert[id_uvw];
						int w = third_vert[id_uvw];
						int id_uv = Arrays.binarySearch(edgeTail, offset[Math.min(u, v)], offset[Math.min(u, v) + 1],
								Math.max(u, v));
						int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u, w)], offset[Math.min(u, w) + 1],
								Math.max(u, w));
						int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(w, v)], offset[Math.min(w, v) + 1],
								Math.max(w, v));

						if ((edge_in_C_k.get(id_uv) == true) && (edge_in_C_k.get(id_uw) == true)
								&& (edge_in_C_k.get(id_vw) == true)) {
							// The triangle is alive...
							// Now, check whether it has enough degree.

							List<Integer> cli = triangle_comneigh.get(id_uvw);
							int cli_len = cli.size();
							int count = 0;
							int j = 0;

							while ((j < cli_len) && (count < k_value)) {
								int z = cli.get(j);
								int id_uz = Arrays.binarySearch(edgeTail, offset[Math.min(u, z)],
										offset[Math.min(u, z) + 1], Math.max(u, z));
								int id_vz = Arrays.binarySearch(edgeTail, offset[Math.min(v, z)],
										offset[Math.min(v, z) + 1], Math.max(v, z));
								int id_wz = Arrays.binarySearch(edgeTail, offset[Math.min(w, z)],
										offset[Math.min(w, z) + 1], Math.max(w, z));
								if ((edge_in_C_k.get(id_uz) == true) && (edge_in_C_k.get(id_vz) == true)
										&& (edge_in_C_k.get(id_wz) == true)) {
									count++;
								}
								j++;
							}
							if (count < k_value) { // we should remove that triangles by removing one of its edges

								update = true;
								member_in_C_k.set(id_uvw, false);

								// Now, we should remove one edge
								double prob1 = edge_prob[id_uv];
								double prob2 = edge_prob[id_uw];
								double prob3 = edge_prob[id_vw];
								int min_index = 0;
								if (prob1 <= prob2) {
									if (prob1 <= prob3) {
										min_index = id_uv;
									} else {
										min_index = id_vw;
									}
								} else {
									if (prob2 <= prob3) {
										min_index = id_uw;
									} else {
										id_uw = id_vw;
									}
								}
								edge_in_C_k.set(min_index, false);
							}
						} else {
							// The triangle should have been removed
							member_in_C_k.set(id_uvw, false);
							update = true;
						}
					}
				}
			}
		}
	}

	
	public double find_density(Set<Integer> edge_list, int num_vertices) {
		double sum = 0.0;
		for (Integer id_e : edge_list) {
			sum += edge_prob[id_e.intValue()];
		}
		
		
		//int num_vertices = vert_list.size();
		return sum / ((double) 0.5 * num_vertices * (num_vertices - 1));
	}

	public void Check_condition(BufferedWriter wg, Set<Integer> triangle_in_Q, BitSet edge_in_Q,int k_value) throws IOException {

		// for each triangle, it stores the number of samples which are deterministic
		// k-nucleus, and contains that triangle
		int[] num_sample_contained_in = new int[num_tria_alive];

		for (int N = 0; N < num_samples; N++) {
			List<Integer> triangle_list_sampleGraph = new ArrayList<Integer>(num_tria_alive);
			boolean nucleus_sample = true;
			// first, project edges (intersection of the edges exist in both Q and the
			// sample graph)
			BitSet edge_in_N = edge_in_sample.get(N);
			BitSet edge_intersection = (BitSet) edge_in_Q.clone();
			// System.out.println(edge_intersection.cardinality());

			// It contains all the edges exist in both N and Q
			edge_intersection.and(edge_in_N);
			// System.out.println(edge_intersection.cardinality());

			for (Integer id_uvw : triangle_in_Q) {

				int u = first_vert[id_uvw];
				int v = second_vert[id_uvw];
				int w = third_vert[id_uvw];

				int id_uv = Arrays.binarySearch(edgeTail, offset[Math.min(u, v)], offset[Math.min(u, v) + 1],
						Math.max(u, v));
				int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u, w)], offset[Math.min(u, w) + 1],
						Math.max(u, w));
				int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(w, v)], offset[Math.min(w, v) + 1],
						Math.max(w, v));

				if ((edge_intersection.get(id_uv) == true) && (edge_intersection.get(id_uw) == true)
						&& (edge_intersection.get(id_vw) == true)) {

					// triangle (u,v,w) is in the sample graph
					triangle_list_sampleGraph.add(id_uvw);
					int count = 0;
					List<Integer> _4clique = triangle_comneigh.get(id_uvw);
					int len = _4clique.size();
					int j = 0;
					while ((j < len) && (count < k_value)) {
						int z = _4clique.get(j);

						int id_uz = Arrays.binarySearch(edgeTail, offset[Math.min(u, z)], offset[Math.min(u, z) + 1],
								Math.max(u, z));
						int id_vz = Arrays.binarySearch(edgeTail, offset[Math.min(v, z)], offset[Math.min(v, z) + 1],
								Math.max(v, z));
						int id_wz = Arrays.binarySearch(edgeTail, offset[Math.min(w, z)], offset[Math.min(w, z) + 1],
								Math.max(w, z));

						if ((edge_intersection.get(id_uz) == true) && (edge_intersection.get(id_vz) == true)
								&& (edge_intersection.get(id_wz) == true)) {
							count++;
							// System.out.println("yes");
						}

						j++;

					}

					if (count < k_value) {
						nucleus_sample = false;
						break;
					}
				}

			}

			if (nucleus_sample == true) {

				for (int j = 0; j < triangle_list_sampleGraph.size(); j++) {
					int id_uvw = triangle_list_sampleGraph.get(j);
					num_sample_contained_in[id_uvw]++;
					// System.out.println(num_sample_contained_in[id_uvw]);

				}
			}
		}

		boolean st = true;
		for (Integer id_uvw : triangle_in_Q) {
			double alpha_k = (double) num_sample_contained_in[id_uvw] / (double) num_samples;
			

			//if (alpha_k < eta2) { // prev 0.0001
			if (alpha_k < eta) { // prev 0.0001
				st = false;
				break;
			}
		}
		if (st == true) {
			// System.out.println("the results are: ");
			Set<Integer> egde_list = new HashSet<Integer>();
			Set<Integer> vertex_list = new HashSet<Integer>();
			
			for (Integer id_uvw : triangle_in_Q) {
				// System.out.print(id_uvw + " ");
				int u = first_vert[id_uvw];
				int v = second_vert[id_uvw];
				int w = third_vert[id_uvw];
				
				int id_uv = Arrays.binarySearch(edgeTail, offset[Math.min(u, v)], offset[Math.min(u, v) + 1],
						Math.max(u, v));
				int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u, w)], offset[Math.min(u, w) + 1],
						Math.max(u, w));
				int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(w, v)], offset[Math.min(w, v) + 1],
						Math.max(w, v));

				egde_list.add(id_uv);
				egde_list.add(id_uw);
				egde_list.add(id_vw);
				
				vertex_list.add(u);
				vertex_list.add(v);
				vertex_list.add(w);
				
				//System.out.print("("+u+","+v+","+w+")" + " ");
				
				//wg.write("("+u+","+v+","+w+")" + " "); 
			}
			
			double density = find_density(egde_list, vertex_list.size());
			double pcc = find_clust_coef(vertex_list, egde_list, triangle_in_Q);
			
			wg.write(density + "\t" + pcc + "\t" +  egde_list.size() + "\t" + vertex_list.size() + "\n");
			//System.out.print(density + " " + pcc);
			
			//for(Integer e: egde_list) {
				//int u = findHead(e);
				//int v = edgeTail[e];
				//wg.write("\""+ author_id[u] + "\"" + "->" + "\""+ author_id[v] + "\"");
				//wg.write("\n");
			//}
			
			//wl.write("\""+ author_id[u] + "\"" + "->" + "\""+ author_id[v] + "\"");
			//wl.write("\n");
			// System.out.println();
			//wg.write("****************************************************************" + "\n");
			 //wg.write("\n");
			 //wg.write("******************" + "\n");
		}
	}
	public double find_clust_coef(Set<Integer> vert_list, Set<Integer> edge_list, Set<Integer> tri_list) {

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
		
		for (Integer id_uvw : tri_list) {
			
			int u = first_vert[id_uvw.intValue()];
			int v = second_vert[id_uvw.intValue()];
			int w = third_vert[id_uvw.intValue()];

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
						if ((support[id_uvz] >= k_val) && (support[id_uwz] >= k_val)
								&& (support[id_vwz] >= k_val)) {
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

	List<Integer> Find_4clique(int u, int v, int w, int[] u_neighbors, int[] v_neighbors, int[] w_neighbors, int u_deg,
			int v_deg, int w_deg) {

		int degree = 0;

		int minlength = 0;
		minlength = Math.min(u_deg, Math.min(v_deg, w_deg));
		List<Integer> clique_list = new ArrayList<Integer>(minlength);

		for (int i = 0, j = 0, k = 0; i < u_deg && j < v_deg && k < w_deg;) {

			if (u_neighbors[i] == v_neighbors[j] && u_neighbors[i] == w_neighbors[k]) { // Find a 4-clique !
//            System.out.println(u + ", " + v + ", " + w  + ", " + u_neighbors[i]); // for checking
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
			
			return findTrinagle_Id(u, v, w);
		}

	}

	public List<Integer> Read_AuthorsOfInterest_AuthorSet_FindId(String filename1, String filename2) throws Exception {

		List<Integer> authors_of_interest_id = new ArrayList<Integer>();
		Set<String> authors_of_interest = new HashSet<String>(); 
		
		BufferedReader br1 = new BufferedReader(new FileReader(filename1));
		String line1 = br1.readLine();

		while (line1 != null) {
			String author = line1.split("\\s+")[0];
			authors_of_interest.add(author);
			line1 = br1.readLine();
		}

		BufferedReader br2 = new BufferedReader(new FileReader(filename2));
		String line2 = br2.readLine();
		
		
		int count =0;
		while (line2 != null) {
			int id = Integer.parseInt(line2.split("\\s+")[0]);
			String author = line2.split("\\s+")[1];
			//author_id[count]= author;
			count++;
			if (authors_of_interest.contains(author)) {
				authors_of_interest_id.add(id);
				
				System.out.println(author + " " + id);
				
			}
			line2 = br2.readLine();
		}
		
		br2.close();
		br1.close();
		return authors_of_interest_id;

	}

	public void ReadFile(String filename) throws Exception {
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

		System.out.println("Number of Triangles with high existence probability: " + t_total_cnt);
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

	public static void main(String[] args) throws Exception {
		//String basename = "DBLP_data_net-proc.w";
		//int precision = 16;
		//double eta = 0.0000000000000001;
		//int k_value = 1;
		//double threshold2 = 0.0000000000000001;
		//String keyword = "data_alg";
		
		
		String basename = args[0];
		int precision = Integer.parseInt(args[1]);
		double eta = Double.parseDouble(args[2]);
		int num_sample = Integer.parseInt(args[3]);
		//double eta2 = Double.parseDouble(args[4]);
		
		//String keyword = args[3];
		//int num_sample = Integer.parseInt(args[4]);
		
		//String basename = "DBLP_data_applic-proc.w";
		
		//String basename = "krogan-proc.w";
		//String basename = "Flickr-proc.w";
		//int precision = 16;
		//double eta = 0.00000000001;
		//double eta = 0.000001;
		//double eta = 0.0000000000001;
		//double eta = 0.1;
		
		//double eta = 0.0000000000000001;
		//int num_sample = 200;
		//double eta2 = 0.1;
		

		Golbal_nuclei_Finding_onSmallGraph_ForExperiments gb = new Golbal_nuclei_Finding_onSmallGraph_ForExperiments(basename, precision, eta,num_sample);
		gb.findConnectedComponnents();
	}

}
