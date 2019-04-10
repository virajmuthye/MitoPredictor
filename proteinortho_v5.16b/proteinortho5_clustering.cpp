/*
 *	Clustering algorithm for Proteinortho
 *	Reads edge list and splits connected components
 *	according to algebraic connectivity threshold
 *
 *	Last updated: 2017/09/20		
 *	Author: Marcus Lechner
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>
#include <stack>
#include <iomanip>
#include <cstdlib>
using namespace std;

// Structs
struct protein {vector<unsigned int> edges; unsigned int species_id; string full_name;};

// Functions
double string2double(string);
void tokenize(const string& , vector<string>& , const string&);
void parse_file(string);
void remove_edge_between(const unsigned int, const unsigned int);
void remove_edge(protein&, const unsigned int);
double getConnectivity(vector<unsigned int>&);
void clear_edges(vector<unsigned int>&);
void partition_graph(void);
void print_group(vector<unsigned int>& , double );
void print_header(void);
void sort_species(void);
vector<unsigned int> get_deg_one (vector<unsigned int>& );

// Parameters
bool param_verbose 		= false;
double param_con_threshold 	= 0.1;
unsigned int debug_level	= 0;
double param_sep_purity 	= 0.75;
string param_rmgraph            = "remove.graph";
unsigned int min_iter 		= 100;		// min number of alg con iterations
unsigned int max_iter 		= 10000;	// max number of alg con iterations

// Globals
unsigned int species_counter = 0;	// Species
unsigned int protein_counter = 0;	// Proteins
vector<string> species;			// Number -> Name
vector<protein> graph;			// Graph containing all protein data
double last_stat = 0;			// For progress stats
unsigned int edges = 0;			// number of edges
ofstream graph_clean;			// File to store graph data
vector<int> reorder_table;		// Tells how proteins/species must be sorted

// TMP Globals
map<string,int> species2id;		// Name -> Number
map<string,int> protein2id;		// Name -> Number


///////////////////////////////////////////////////////////
// Main
///////////////////////////////////////////////////////////
void printHelp() {
	cerr << "Proteinortho5 - Spectral partitioning algorithm" << endl;
	cerr << "-----------------------------------------------" << endl;
	cerr << "This tool is part of Proteinortho" << endl;
	cerr << "" << endl;
	cerr << "Usage:   proteinortho5_clustering [OPTIONS] graph_files..." << endl;
	cerr << "Options: -verbose        report progress" << endl;
	cerr << "         -conn double     threshold for connectivity [0.1]" << endl;
	cerr << "         -rmgraph STRING output file for graph" << endl;
}

int main(int argc, char *argv[]) {
	if (argc <= 1) {
		printHelp();
		return EXIT_FAILURE;
	}

	try {
	// Read parameters
	int paras;
	vector<string> files;
	for (paras = 1; paras < argc; paras++) {
		string parameter = string(argv[paras]);
		if (parameter.substr(0, 1) != "-") {
			files.push_back(parameter);
		}
		else if (parameter == "-verbose") {
			paras++;
			if (string2double(string(argv[paras]))) {
				param_verbose = true;
			}
		}
		else if (parameter == "-conn") {
			paras++;
			param_con_threshold = string2double(string(argv[paras]));
		}
		else if (parameter == "-purity") {
			paras++;
			param_sep_purity = string2double(string(argv[paras]));
		}
		else if (parameter == "-debug") {
			paras++;
			debug_level = int(string2double(string(argv[paras])));
		}
		else if (parameter == "-rmgraph") {
			paras++;
			param_rmgraph = string(argv[paras]);
		}
		else {
			printHelp();
			cerr << endl << "Sorry, unknown option '" << string(argv[paras]) << "'!" << endl;
			return EXIT_FAILURE;
		}
	}
	
	if (debug_level > 0) cerr << "[DEBUG] Debug level " << debug_level << endl;

	// Parse files
	for (vector<string>::iterator it=files.begin() ; it != files.end(); it++) {
		if (debug_level > 0) cerr << "[DEBUG] Parsing file " << *it << endl;
		parse_file(*it);
		if (debug_level > 0) cerr << "[DEBUG] I know " << species_counter <<  " species with " << protein_counter << " proteins and " << edges << " edges in sum" << endl;
	};

	// Free memory
	files.clear();
	species2id.clear();
	protein2id.clear();

	// Stats
	if (param_verbose) cerr << species_counter << " species" << endl << protein_counter << " paired proteins" << endl << edges << " bidirectional edges" << endl;

	// Prepare sort of output
	if (debug_level > 0) cerr << "[DEBUG] Sorting known species" << endl;
	sort_species();

	// Write output header
	print_header();							

	// Open graph-removal file
	graph_clean.open(param_rmgraph.c_str());

	// Clustering
	if (debug_level > 0) cerr << "[DEBUG] Clustering" << endl;
	partition_graph();
	graph_clean.close();
	}
	catch( string& error ) {
	      cout << "Error: " << error << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////
// Debug
///////////////////////////////////////////////////////////
void print_graph() {
	// Print graph
	for (unsigned int i = 0; i < graph.size(); i++) {
		for (unsigned int j = 0; j < graph[i].edges.size(); j++) {
//			cerr << i << " -> " << graph[i].edges[j] << " \t \t ";
			cerr << graph[i].full_name << " -> " << graph[graph[i].edges[j]].full_name << endl;
		}
	}
}

void print_edgelist(protein& protein) {
	for (vector<unsigned int>::iterator ita = protein.edges.begin(); ita != protein.edges.end(); ita++) {
		cerr << graph[*ita].full_name << " ";
	}
	cerr << endl;
}

///////////////////////////////////////////////////////////
// Output
///////////////////////////////////////////////////////////
// Sort
void sort_species(void) {
	reorder_table.reserve(species_counter);
	vector<string> species_sorted (species_counter);
	copy(species.begin(), species.end(), species_sorted.begin());
	sort(species_sorted.begin(), species_sorted.end());
	// find new locations (not efficent but list is small)	
	for (unsigned int i = 0; i < species_counter; i++) {
		for (unsigned int j = 0; j < species_counter; j++) {
			if (species[i] == species_sorted[j]) {reorder_table[j] = i; continue;}
		}
	}
}


// Progress stats
void stats(double i, double size) {
	if (!param_verbose) return;
	double stat = double(i/size*100);
	if (last_stat * 1.01 < stat) {
		last_stat = stat;
		cerr << "\r" << "                          ";
		cerr << "\r" << "Clustering: " << setprecision (2) << fixed << stat << "%";
//		cerr << "Clustering: " << stat << "%" << endl;
	}
}

// Header with species names
void print_header() {
	cout << "# Species\tGenes\tAlg.-Conn.";
	for (unsigned int i = 0; i < species_counter; i++) {
		cout << "\t" << species[reorder_table[i]];
	}
	cout << endl;
}

// Group formatting
void print_group(vector<unsigned int>& nodes, double connectivity) {
	vector<string> line(species_counter,"*");	// Output vector

	unsigned int species_number = 0;
	// For each protein in group
	for (unsigned int i = 0; i < nodes.size(); i++) {
		unsigned int current_protein = nodes[i];
		unsigned int current_species = graph[current_protein].species_id;
		if (line[current_species].compare("*") == 0) {
			line[current_species] = graph[current_protein].full_name;
			species_number++;
		}
		else {
			line[current_species].append(","+graph[current_protein].full_name);
		}
	}

	cout << species_number << "\t" << nodes.size() << "\t" << setprecision (3) << connectivity;

	// List group data
	for (unsigned int i = 0; i < species_counter; i++) {
		
		// sort line if multi protein
		vector<string> fields;
		tokenize(line[reorder_table[i]], fields, ",");
		if (fields.size() > 1) {
			sort(fields.begin(), fields.end());
			line[reorder_table[i]] = fields[0];
			for (unsigned int k = 1; k < fields.size(); k++) {
				line[reorder_table[i]].append(","+fields[k]);
			}
		}

		// output
		cout << "\t" << line[reorder_table[i]];
	}
	cout << endl;
}


///////////////////////////////////////////////////////////
// Major partioning algorithm
///////////////////////////////////////////////////////////
void partition_graph() {
	vector<bool> done = vector<bool> (protein_counter, false);	// Keep track on what was done
	srand(12345);							// Init random number generator

	// For each protein (no increment here, we might redo protein i)
	for (unsigned int protein_id = 0; protein_id < graph.size();) {
//		cerr << "Protein ID: " << protein_id << endl;
		// We were here already
		if (done[protein_id])		{protein_id++;	continue;}

		// Mark protein as done
		done[protein_id] = true;

		// Init todo list with this protein
		stack<unsigned int> todo;				// Coloring stack
		todo.push(protein_id);

		// Collect members of the current connected component
		vector<unsigned int> current_group;			// Keep track of group
		while (!todo.empty()) {
			// Get next protein & remove it from stack
			unsigned int next = todo.top();
			todo.pop();
			// Add protein to current group
			current_group.push_back(next);

			// Add targets to todo list
			for (unsigned int i = 0; i < graph[next].edges.size(); i++) {
				unsigned int target = graph[next].edges[i];
				// Mark & add if unknown yet
				if (!done[target]) {
					done[target] = true;
					todo.push(target);
				}
			} // For each target
		} // While todo list is not empty		

		// Do not report singles
		if (current_group.size() < 2) {
			if (debug_level > 1) cerr << "[DEBUG] Single gene skipped" << endl;
			protein_id++;
			continue;
		}

		// Connectivity analysis
		if (debug_level > 0) cerr << "[DEBUG] Calculating connectivity of a group. " << current_group.size() << " proteins (ID: " << protein_id << ")" << endl;
		if (current_group.size() > 16000) cerr << "[WARN] Current connected component contains " << current_group.size() << " proteins. That might be way to much for me. Try raising the e-value!" << endl;
		double connectivity = getConnectivity(current_group);
		if (debug_level > 0) cerr << "[DEBUG] Connectivity was " << connectivity << endl;

		if (connectivity < param_con_threshold && current_group.size() > 3) {
			if (debug_level > 0) cerr << "[DEBUG] Reiterating" << endl;
			// Split groups is done in getConnectivity function
			// Reset flags and repeat without incrementing protein counter
			for (unsigned int i = 0; i < current_group.size(); i++) done[current_group[i]] = false;
			continue;
		}
	
		// Output
		if (connectivity >= param_con_threshold) {print_group(current_group,connectivity);}

		// Print stats
		stats(protein_id,protein_counter);
		
		// Clean up and go on with the next connected component
		clear_edges(current_group);
		protein_id++;
	} // For each protein (mainloop)
	stats(1,1);
	if (param_verbose) cerr << "\r" << "Done                       " << endl;
}

///////////////////////////////////////////////////////////
// Algebraic connectivity functions
///////////////////////////////////////////////////////////
// Return maximum degree of given protein_ids
unsigned int max_deg(vector<unsigned int>& nodes) {
	unsigned int max = 0;
	for (unsigned int i = 0; i < nodes.size(); i++) {
		unsigned int degree = graph[nodes[i]].edges.size();
		if (degree > max) max = degree;
	}
	return max;
}

// Return nodes with degree 1 of given protein_ids
vector<unsigned int> get_deg_one (vector<unsigned int>& nodes) {
	vector<unsigned int> one;
	for (unsigned int i = 0; i < nodes.size(); i++) {
		unsigned int degree = graph[nodes[i]].edges.size();
		if (degree == 1) {one.push_back(nodes[i]);}
	}
	return one;
}

// Generate random vector x of size size
vector<double> generate_random_vector(const unsigned int size) {
	vector<double> x(size);
	for (unsigned int i = 0; i < size; i++) {
	  x[i] = (double)(rand() % 999+1)/1000;	// 1 bis 99
		// at least one value must be different from the others but still within 0 and 1
		if (i > 0 && x[i] == x[i-1]) {
			x[i] /= 3;
		}
//		cerr << x[i] << endl;
	}
	return x;
}

// Generate random vector x of size size
vector<double> generate_random_vector_old(const unsigned int size) {
	vector<double> x(size);
	for (unsigned int i = 0; i < size; i++) {
	  x[i] = (double)rand()/RAND_MAX;
		// at least one value must be different from the others but still within 0 and 1
		if (i > 0 && x[i] == x[i-1]) {
			x[i] += 0.1;
			if (x[i] > 1) {
				x[i] -= 0.2;
			}
		}
//		cerr << x[i] << endl;
	}
	return x;
}

// determine new X, Formula (1)
vector<double> get_new_x(vector<double> x, vector<unsigned int>& nodes, map<unsigned int,unsigned int>& mapping) {
	vector<double> x_new(x.size());
	// go through all nodes (rows of A)
	for (unsigned int i = 0; i < nodes.size(); i++) {
		unsigned int node = nodes[i];	// node x is node y here
		x_new[i] = 0;
		// go through adjacency list of node (cols of A_i)
		for (unsigned int j = 0; j < graph[node].edges.size(); j++) {
			// y points to z, so take entry z from x
			unsigned int abs_target = graph[node].edges[j];
			unsigned int rel_target = mapping[abs_target];
			// -> Sum over x entries
			x_new[i] += x[rel_target];
		}
	}
	return x_new;
}

// Make vector x orthogonal to 1, Formula (2)
vector<double> makeOrthogonal(vector<double> x) {
	double sum = 0;
	for (unsigned int i = 0; i < x.size(); i++) {sum += x[i];}
	double average = sum/x.size();
	for (unsigned int i = 0; i < x.size(); i++) {x[i] -= average;}
	return x;
}

// Normalize vector x, Formula (4)
vector<double> nomalize(vector<double> x, double *length) {
	double sum = 0;
	for (unsigned int i = 0; i < x.size(); i++) {sum += x[i]*x[i];}
	*length = sqrt(sum);
	if (*length == 0) {*length = 0.000000001;}	// ATTENTION not 0!
	for (unsigned int i = 0; i < x.size(); i++) {x[i] /= *length;}
	return x;
}

// Qx, Formula (5)
vector<double> getY(double max_degree, vector<double> x_hat, vector<double> x_new, vector<unsigned int>& nodes){
	// (2*maxdeg - grad_node_i ) * x_hat_i + new_x_i
	for (unsigned int i = 0; i < x_hat.size(); i++) {
		x_hat[i] *= (2*max_degree - graph[nodes[i]].edges.size());
		x_hat[i] += x_new[i];
	}
	return x_hat;
}


// Remove edges connectiong two groups a and b
void removeExternalEdges(map<unsigned int,bool>& a) {
//		cerr << "+#" << endl;
//		for (map<unsigned int,bool>::iterator it = a.begin(); it != a.end(); it++) {
//			unsigned int protein = it->first;
//			cerr << protein << endl;
//		}
//		cerr << "#-" << endl;
		
		// For each protein in a
		for (map<unsigned int,bool>::iterator it = a.begin(); it != a.end(); it++) {
			unsigned int protein = it->first;
			// For each target
			vector<unsigned int> cleaned_edges;
			bool swap = false;
			for (vector<unsigned int>::iterator ita = graph[protein].edges.begin(); ita != graph[protein].edges.end(); ita++) {
				// If it is not present the own group, set flag
				if (a.find(*ita) == a.end()) {
					// 5.16 store name AND species in removal list
					graph_clean << graph[protein].full_name << "\t" << species[graph[protein].species_id] << "\t" << graph[*ita ].full_name << "\t" << species[graph[*ita ].species_id] << endl; // Improved graph cleaning
					swap = true;
				}
				// Otherwise, add it to the new edge list
				else {
					cleaned_edges.push_back(*ita);
				}
			}
			// If changes were made, swap edge list with new one		
			if (swap) {
				cleaned_edges.swap(graph[protein].edges);
			}
		}
}


// Remove edges connectiong two groups a and b
void removeExternalEdges_old(map<unsigned int,bool>& a) {
		// For each protein in a
		for (map<unsigned int,bool>::iterator it = a.begin(); it != a.end(); it++) {
			unsigned int protein = it->first;
			// For each target
			vector<unsigned int> cleaned_edges;
			bool swap = false;
			for (vector<unsigned int>::iterator ita = graph[protein].edges.begin(); ita != graph[protein].edges.end(); ita++) {
				// If it is not present the own group, set flag
				if (a.find(*ita) == a.end()) {
					graph_clean << graph[protein].full_name << "\t" << graph[*ita ].full_name << endl; // Improved graph cleaning
					swap = true;
				}
				// Otherwise, add it to the new edge list
				else {
					cleaned_edges.push_back(*ita);
				}
			}
			// If changes were made, swap edge list with new one		
			if (swap) {
				cleaned_edges.swap(graph[protein].edges);
			}
		}
}

// Split connected component according to eigenvector
void splitGroups(vector<double>& y, vector<unsigned int>& nodes){

//	cerr << "Mission to split groups:" << endl << "####################" << endl;		///

//	cerr << "Vector (" << nodes.size() << ")" << endl;
//	for (unsigned int i = 0; i < nodes.size(); i++) {
//		cerr << graph[nodes[i]].full_name << "\t" << nodes[i] << "\t" << y[i] << endl;
//	}


	// Remove tree like structures in the beginning
//	vector<unsigned int> one = get_deg_one(nodes);
//	if (one.size() > 0) {
//		cerr << "Tree " << endl;		///
//		map<unsigned int,bool> tree;
//		for (unsigned int i = 0; i < one.size(); i++) {
//			{tree[one[i]] = true;}
//		}
//		removeExternalEdges(tree);
//		return;
//	}

	// Store data about two groups (Zero cannot be assigned with certainty)
	map<unsigned int,bool> groupA, groupB, groupZero;
	for (unsigned int i = 0; i < y.size(); i++) {
		unsigned int style = 1;
		if (y[i] < 0) style = 2;
		if (abs(y[i]) < param_sep_purity) style = 0;
		if (debug_level > 2) cerr << "[DEBUG L3] " << i << " = " << y[i] << " -> " << style << endl;
		if      (abs(y[i]) < param_sep_purity)	{groupZero[nodes[i]] = true;} // cerr << graph[nodes[i]].full_name << " {color:#95cde5}" << endl; }
		else if (y[i] < 0) 				{groupA[nodes[i]] = true;} // cerr << graph[nodes[i]].full_name << " {color:#95cde5}" << endl; }
		else              				{groupB[nodes[i]] = true;} // cerr << graph[nodes[i]].full_name << " {color:#b01700}" << endl; }
	}

//	cerr << groupA.size() << " -- " << groupB.size() << endl;

	// Catch error in laplacien calcs
	if ((groupA.size() == 0 || groupB.size() == 0) && groupZero.size() == 0) {
		throw string("Failed to partition subgraph! This might lead to an infinit loop. Please submit the .blastgraph file to lechner@staff.uni-marburg.de to help fixing this issue.");
	}

	removeExternalEdges(groupZero);
	removeExternalEdges(groupA);
	removeExternalEdges(groupB);
}

double getConnectivity(vector<unsigned int>& nodes) {
	// Special case quicky
	if (nodes.size() == 2) {
		if (debug_level > 1) cerr << "[DEBUG L2] Shortcut for two-piece component -> conn = 1" << endl;
		return 1;
	}
	// quicky end

	// Get max degree of nodes
	unsigned int max_degree = max_deg(nodes);
	
	// Compress value range / map data
	map<unsigned int,unsigned int> mapping;
	for (unsigned int i = 0; i < nodes.size(); i++) {mapping[nodes[i]] = i;}
	
	// Init randomized x 
	vector<double> x = generate_random_vector(nodes.size());

	// Orthogonalize + normalize vector + get initial lenght
	double last_length, current_length = 0;
	vector<double> x_hat = makeOrthogonal(x);
	vector<double> norm = nomalize(x_hat, &last_length);

	// Repeat until difference < epsilon
	unsigned int iter = 0;							// catch huge clustering issues by keeping track here
	
	while(iter++ < max_iter) { 
		if (debug_level > 2 && iter%100 == 0) cerr << "[DEBUG L2] Step " << iter << " / " << max_iter << endl;
		last_length = current_length;
		// Get a new x
		x = get_new_x(norm, nodes, mapping);
		// Get y
		vector<double> y = getY(max_degree,norm,x,nodes);
		// Orthogonalize
		x_hat = makeOrthogonal(y);
		// Get lenght (lambda) & normalize vector
		norm = nomalize(x_hat, &current_length);
///		cerr << "I " << iter << ": " << abs(current_length-last_length) << endl; 
		if (abs(current_length-last_length) < 0.0001 && iter >= min_iter) break;	// min 100 iterations (prevent convergence by chance), converge to 1e-6
	}
///		cerr << nodes.size() << " nodes done after " << iter << " iterations" << endl;
	if (debug_level > 0) cerr << "[DEBUG] " << iter << " / " << max_iter << " iterations required (error is " << abs(current_length-last_length) << ")" << endl;

	double connectivity = (-current_length+2*max_degree)/(nodes.size());

//	cerr << nodes.size() << " " << connectivity << endl;

///	if (connectivity > 1) {cerr << "DIE" << endl; throw "XX";}

	// Split groups if connectivity is too low, remove tree like structures that might have arosen
	if (connectivity < param_con_threshold) {
///		cerr << "Conn: " << connectivity << endl; ///
		splitGroups(x_hat, nodes);
	}
	
	return connectivity;
}

///////////////////////////////////////////////////////////
// Basic Graph functions
///////////////////////////////////////////////////////////
// Remove all edges from the given list of protein ids
void clear_edges(vector<unsigned int>& nodes) {
	for (unsigned int i = 0; i < nodes.size(); i++) graph[nodes[i]].edges.clear();
}

// Remove edge between two
void remove_edge_between(const unsigned int a_id, const unsigned int b_id) {
	remove_edge(graph[a_id], b_id);
	remove_edge(graph[b_id], a_id);
}

// Remove edge from protein
void remove_edge(protein& node, const unsigned int remove_id) {
	// Search for element in edge list
	vector<unsigned int>::iterator element = find(node.edges.begin(), node.edges.end(), remove_id);
	// Not found (something is wrong!)
	if (element == node.edges.end()) throw string("Element could not be found in edge list");
//	// Remove element
	node.edges.erase(element);
}

///////////////////////////////////////////////////////////
// File parser
///////////////////////////////////////////////////////////
void parse_file(string file) {
	if (param_verbose) cerr << "Reading " << file << endl;
	string line;
	ifstream graph_file(file.c_str());
	if (graph_file.is_open()) {
		// For each line
		string file_a = "";	unsigned int file_a_id = 0;
		string file_b = "";	unsigned int file_b_id = 0;
		while (!graph_file.eof()) {
			getline(graph_file, line);
			vector<string> fields;
			tokenize(line, fields, "\t");
			// Header line
			if (fields.size() == 2 && fields[0].substr(0, 1) == "#") {
				file_a = fields[0].substr(2, fields[0].size()-2);
				file_b = fields[1];

				if (file_a == "file_a" && file_b == "file_b") continue;	// Init Header

				// Map species a if not know so far
				if (species2id.find(file_a) == species2id.end())	{
						species.push_back(file_a);
//						cerr << "Species " << species_counter << ": " << file_a << endl;
						species2id[file_a] = species_counter++;
				}
				// Map species b if not know so far
				if (species2id.find(file_b) == species2id.end())	{
//						cerr << "Species " << species_counter << ": " << file_b << endl;
						species.push_back(file_b);
						species2id[file_b] = species_counter++;
				}

				file_a_id = species2id[file_a];
				file_b_id = species2id[file_b];
			}
			// Data line
			else if ((fields.size() == 6 || fields.size() == 8) && fields[0].substr(0, 1) != "#") {
				// a b e1 b1 e2 b2 score

//				cerr << protein_counter << ": " << fields[0] << " <-> " << fields[1] << endl;

				// 5.16 deal with duplicated IDs by adding file ID to protein ID
				string ida = fields[0];
				string idb = fields[1];
				fields[0] += " "; fields[0] += file_a_id;
				fields[1] += " "; fields[1] += file_b_id;

				// 5.16 do not point to yourself
				if (!fields[0].compare(fields[1])) {continue;}

				// A new protein
				if (protein2id.find(fields[0]) == protein2id.end())	{
					protein a;
					a.full_name	= ida;
					a.species_id	= file_a_id;
					protein2id[fields[0]] = protein_counter++;
					graph.push_back(a);
				}
				if (protein2id.find(fields[1]) == protein2id.end())	{
					protein b;
					b.full_name	= idb;
					b.species_id	= file_b_id;
					protein2id[fields[1]] = protein_counter++;
					graph.push_back(b);
				}

				// Add link to graph (reciprocal)					
				unsigned int a_id = protein2id[fields[0]];
				unsigned int b_id = protein2id[fields[1]];

				graph[a_id].edges.push_back(b_id);
				graph[b_id].edges.push_back(a_id);
				edges++;
			}
		}
	}
	else {
		throw string("Could not open file ") + file;
	}
}



///////////////////////////////////////////////////////////
// Misc functions
///////////////////////////////////////////////////////////
// Convert string to double
double string2double(string str) {
	istringstream buffer(str);
	double value;
	buffer >> value;
	return value;
}

// Split a string at a certain delim
void tokenize(const string& str, vector<string>& tokens,
		const string& delimiters = "\t") {
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos) {
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}
