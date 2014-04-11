//=======================================================================
// Path Ranker 
// 
// Author: Ichigaku Takigawa
//   a2ps --no-header --landscape --columns=1 --font-size=3.6
//   output.log -o test.ps
//=======================================================================

#include <deque>
#include <utility>
#include <fstream>
#include <iterator>
#include <math.h>


// boost 1.33.1 required
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/pending/relaxed_heap.hpp>
#include <boost/pending/indirect_cmp.hpp>

#define R_NO_REMAP

#include "init.h"
#include "hme3m.h"
using namespace boost;

//  ****** base graph ******  //

typedef adjacency_list<
  vecS,          // for out-edges
  vecS,           // for vertices
  bidirectionalS, // bidirectional graph
  // vertex properties
  property<vertex_name_t, string>,
  // edge properties
  property<edge_weight_t, double,
       property<edge_name_t, string>
       >
  > Graph;

//  ****** type definitions ******  //

typedef graph_traits<Graph>::vertices_size_type   SizeType;
typedef graph_traits<Graph>::vertex_descriptor    Vertex;
typedef graph_traits<Graph>::edge_descriptor      Edge;
typedef pair<Vertex,Vertex>                       VertexPair;
typedef graph_traits<Graph>::vertex_iterator      VertexIter;
typedef graph_traits<Graph>::edge_iterator        EdgeIter;
typedef graph_traits<Graph>::in_edge_iterator     InEdgeIter;
typedef graph_traits<Graph>::adjacency_iterator   AdjacencyIter;

typedef property_map<Graph, vertex_name_t>::type  VertexNameMap;
typedef property_map<Graph, edge_name_t>::type    EdgeNameMap;
typedef property_map<Graph, edge_weight_t>::type  EdgeWeightMap;
typedef property_map<Graph, vertex_index_t>::type VertexIndexMap;

typedef double   WeightType;
typedef unsigned LengthType;

Graph& get_graph(pair<Graph,VertexPair>& _data){
  return _data.first;
}

Vertex& get_start_vertex(pair<Graph,VertexPair>& _data){
  return _data.second.first;
}

Vertex& get_end_vertex(pair<Graph,VertexPair>& _data){
  return _data.second.second;
}

ifstream& operator>>(ifstream& ifs, VertexPair& vpair){
    ifs >> vpair.first >> vpair.second;
    return ifs;
}

template <typename IndexedSequence>
bool equal_first_i_nodes(const IndexedSequence& seq1,
             const IndexedSequence& seq2,
             const int i) {
  // seq1[0]...seq1[i] and seq2[0]...seq2[i] equal?
  for(int j=0; j<=i; ++j){
    if(seq1[j] != seq2[j]){ return false; }
  }
  return true;
}

struct st_path_with_deviation {
  deque<Vertex> sequence;
  WeightType    score;
  Vertex        deviation;
};

inline bool compare(const struct st_path_with_deviation& x, 
            const struct st_path_with_deviation& y){
  return x.score < y.score;
}

//option
extern char *optarg;
extern int  optind, opterr;

//  ****** dijkstra algorithm ******  //
template<typename PredecessorMap, typename DistanceMap, typename Compare>
void dijkstra_algorithm(Graph& g, const Vertex& s,
            PredecessorMap predecessor,
            DistanceMap distance,
            Compare compare){
  Vertex u, v;
  SizeType n = num_vertices(g);
  EdgeWeightMap  weight = get(edge_weight, g);

  typedef indirect_cmp<DistanceMap, Compare> IndirectCmp;
  IndirectCmp icmp(distance, compare);

  relaxed_heap<Vertex, IndirectCmp> minheap(n,icmp);
  
  VertexIter vitr,vend;
  for(tie(vitr,vend)=vertices(g); vitr != vend; ++vitr) {
    minheap.push(*vitr);
  }

  WeightType rhs;
  AdjacencyIter aitr, aend;

  while(not minheap.empty()){
    // extract-min
    u = minheap.top();
    minheap.pop();

    for(tie(aitr,aend)=adjacent_vertices(u,g); aitr != aend; ++aitr){
      v   = *aitr;

      rhs = distance[u] + weight[edge(u,v,g).first];

      // edge relaxation
      if(distance[v] > rhs){
            distance[v] = rhs;
            minheap.update(v);
            predecessor[v] = u;
       }
    }
  }
}

struct st_path_with_deviation
st_shortest_path(const Vertex s, const Vertex t, Graph& g){
  SizeType n_vertices = num_vertices(g);
  vector<Vertex>     predecessor(n_vertices);
  vector<WeightType> distance(n_vertices, numeric_limits<WeightType>::max());
  vector<LengthType> pathlength(n_vertices, numeric_limits<int>::max());

  distance[s] = 0;
  dijkstra_algorithm(g, s, &predecessor[0], &distance[0], std::less<WeightType>());

  struct st_path_with_deviation p;
  if(distance[t] != numeric_limits<WeightType>::max()){
    for(Vertex v = t; v != s; v = predecessor[v]){
      p.sequence.push_front(v);
    }
    p.sequence.push_front(s);
  }
  
  p.score = distance[t];
  p.deviation = 0;
  
  return p;
}

#undef R_NO_REMAP

pair<Graph,VertexPair> R_get_st_graph_from(SEXP nodes, SEXP edges,SEXP edge_weights);

SEXP store_path_R(deque<Vertex>& st_path, Graph& g, double p_score);

SEXP pathranker(SEXP node_list, SEXP edge_list, SEXP edge_weights, SEXP rk,SEXP minpathsize)
{
  Vertex s,t;
  Graph g;
  unsigned K =(unsigned)(REAL(rk)[0]);
  unsigned mps = (int)(REAL(minpathsize)[0]);
  SEXP new_path,all_paths;
    
  Rf_protect(all_paths = Rf_allocVector(VECSXP,K));

  pair<Graph,VertexPair> data = R_get_st_graph_from(node_list, edge_list, edge_weights);
    
  g = get_graph(data);
  s = get_start_vertex(data);
  t = get_end_vertex(data);
    
  if(s == numeric_limits<Vertex>::max() or t == numeric_limits<Vertex>::max()){
    Rprintf("No vertex start or end vertices found.");
    return R_NilValue;
  }
    
  // ======= Yen-Lawler algorithm ==================================
  // [1] Jin Y. Yen,
  //     Finding the k shortest loopless paths in a network,
  //     Management science, Vol.17, No. 11, July, 1971. 712-716
  // [2] E. L. Lawler,
  //     A procedure for computing the k best solutions to discrete
  //     optimization problems and its application
  //     to the shortest path problem,
  //     Management science, vol. 18, no. 7, 1972. 401-405
  // ===============================================================
 
  deque<struct st_path_with_deviation> result, candidates;
    
  Graph g_copy;
  deque<Vertex> prefix;

  WeightType    accumulated_score;
  EdgeWeightMap edge_weight_map = get(edge_weight, g);
    
  struct st_path_with_deviation p, q;

  candidates.push_back(st_shortest_path(s, t, g));
  
  unsigned k = 0;
  unsigned rsize = 0;
  while (rsize < K) {
    // check if paths exist
    if(candidates.empty()){
      //cout << "No more paths." << endl;
      break;
    }
    
    // output k-th shortest loopless path
    sort(candidates.begin(), candidates.end(), compare);
    if(candidates.size() > K-rsize+1){
      candidates.resize(K-rsize+1);
    }
    
    p = candidates.front();
    candidates.pop_front();
    result.push_back(p);
    
    if(p.score==numeric_limits<WeightType>::max()){
      //cout << "No more paths." << endl;
      break;
    }

    // I want to make sure the path involves more than 1 gene
    // so the path distance must be greater than 2 times the s-> gene and gene->t
    // edges have the same weight.
    Edge e_temp = edge(p.sequence[0],p.sequence[1],g).first;
    double fw = get(edge_weight, g)[e_temp];

    if (p.sequence.size() > mps && p.score > 2*fw) {
      new_path = store_path_R(p.sequence,g,p.score);
      SET_VECTOR_ELT(all_paths,rsize,new_path);
      rsize = rsize + 1;
    }
    k = k + 1;

    prefix.clear();
    accumulated_score = 0;
    
    for(unsigned i=0; i < p.sequence.size()-1; ++i){
      if(i >= p.deviation){
          // copy of the graph for editing
          g_copy = g;
    
          // edge deletion
          for(unsigned j=0; j < result.size(); ++j){
            if(equal_first_i_nodes(result[j].sequence, p.sequence ,i)){
              remove_edge(result[j].sequence[i], result[j].sequence[i+1], g_copy);
            }
          }             
          // vertex deletion
          for(unsigned j=0; j < i; ++j){
            clear_vertex(p.sequence[j], g_copy);
          }
    
          // compute shortest path
          q = st_shortest_path(p.sequence[i], t, g_copy);

          // concatenate the prefix to the current path
          if(q.score!=numeric_limits<WeightType>::max()){
            if(not q.sequence.empty()){
              q.sequence.insert(q.sequence.begin(), prefix.begin(), prefix.end());
              q.score += accumulated_score;
            }
            q.deviation = i;
            candidates.push_back(q);
          }
        }
      
        prefix.push_back(p.sequence[i]);
        accumulated_score += edge_weight_map[edge(p.sequence[i], p.sequence[i+1],g).first];
    }
  }
   
  Rf_unprotect(1);
  return all_paths;
}

pair<Graph,VertexPair> R_get_st_graph_from(SEXP nodes, SEXP edges,SEXP edge_weights) {

  Vertex s,t;
  s = numeric_limits<Vertex>::max();
  t = numeric_limits<Vertex>::max();
  // read the vertex file
  int counter=0;
  string line;
  vector<string> name;              // names of vertices
  SEXP from_idxs,to_idxs,labels;

  for (int i = 0;i < LENGTH(nodes);i = i + 1) {
    line = CHAR(STRING_ELT(nodes,i));
    name.push_back(line);
    if(line=="s"){
      s = counter;
    }else if(line=="t"){
      t = counter;
    }
    counter++;
  }
  
  SizeType n_vertices = name.size(); // # of vertices
   
  // create a graph object
  Graph g(n_vertices);
  
  // resister vertex names
  VertexIter v_iter, v_end;
  VertexNameMap vertex_name_map     = get(vertex_name, g);
  
  for ( tie(v_iter,v_end)=vertices(g); v_iter != v_end; ++v_iter ){
    vertex_name_map[*v_iter] = name[*v_iter];
  }
  
  // read the edgelist file
  VertexPair  v_pair;
  Edge        e_temp;
  //WeightType  d_value;
  string      edge_label;
  
  EdgeNameMap   edge_name_map   = get(edge_name, g);
  
  // define R types
  from_idxs = VECTOR_ELT(edges,0);
  to_idxs = VECTOR_ELT(edges,1);
  labels = VECTOR_ELT(edges,2);

  for (int i = 0;i < LENGTH(from_idxs);i = i + 1) {
    v_pair.first = (int)(REAL(from_idxs)[i] - 1);
    v_pair.second = (int)(REAL(to_idxs)[i] - 1);

    add_edge(v_pair.first, v_pair.second , REAL(edge_weights)[i] , g);
    e_temp = edge(v_pair.first, v_pair.second , g).first;
    edge_name_map[e_temp] = CHAR(STRING_ELT(labels,i));
  }
  
  pair<Graph,VertexPair> ret;
  ret.first = g;
  ret.second.first  = s;
  ret.second.second = t;
  
  return ret;
}

SEXP store_path_R(deque<Vertex>& st_path, Graph& g, double p_score) {
  SEXP new_path,dimnames;
  SEXP genes,compounds,weights,distance;  

  st_path.pop_back();	st_path.pop_front();  //delete "s" & "t" nodes


  if(st_path.empty()){
    return R_NilValue;
  }
  VertexNameMap vertex_name_map = get(vertex_name, g);
  EdgeNameMap   edge_name_map   = get(edge_name, g);
  EdgeWeightMap edge_weight_map = get(edge_weight, g);
  Edge e_temp;
     
  Rf_protect(genes = Rf_allocVector(STRSXP,st_path.size()));
  Rf_protect(compounds = Rf_allocVector(STRSXP,st_path.size() - 1));
  Rf_protect(weights = Rf_allocVector(REALSXP,st_path.size()-1));
  Rf_protect(distance = Rf_allocVector(REALSXP,1));
  REAL(distance)[0] = p_score;

  SET_STRING_ELT(genes,0,Rf_mkChar(vertex_name_map[st_path[0]].c_str()));
  for (unsigned i = 0; i < st_path.size()-1; ++i) {
    e_temp = edge(st_path[i],st_path[i+1],g).first;
    
    SET_STRING_ELT(compounds,i,Rf_mkChar(edge_name_map[e_temp].c_str()));
    SET_STRING_ELT(genes,i+1,Rf_mkChar(vertex_name_map[st_path[i+1]].c_str()));
    REAL(weights)[i] = edge_weight_map[e_temp];
  }

  Rf_protect(new_path = Rf_allocVector(VECSXP, 4));
  Rf_protect(dimnames = Rf_allocVector(VECSXP,4));

  SET_VECTOR_ELT(new_path,0,genes);
  SET_VECTOR_ELT(dimnames,0,Rf_mkString("genes"));

  SET_VECTOR_ELT(new_path,1,compounds);
  SET_VECTOR_ELT(dimnames,1,Rf_mkString("compounds"));
  
  SET_VECTOR_ELT(new_path,2,weights);
  SET_VECTOR_ELT(dimnames,2,Rf_mkString("weights"));

  SET_VECTOR_ELT(new_path,3,distance);
  SET_VECTOR_ELT(dimnames,3,Rf_mkString("distance"));

  Rf_setAttrib(new_path,R_NamesSymbol,dimnames);

  Rf_unprotect(6);

  return new_path;
}


