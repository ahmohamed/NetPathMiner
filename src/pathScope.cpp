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

//#define MATHLIB_STANDALONE
#include <math.h>
#include <float.h>

//init
#include "init.h"

// Disable warnings in boost header code.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-local-typedef"
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif


// boost 1.33.1 required
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/pending/relaxed_heap.hpp>
#include <boost/pending/indirect_cmp.hpp>

#pragma GCC diagnostic pop

#ifdef __clang__
#pragma clang diagnostic pop
#endif


#define R_NO_REMAP

using namespace boost;

//  ****** base graph ******  //

typedef adjacency_list<
  vecS,          // for out-edges
  vecS,           // for vertices
  bidirectionalS, // bidirectional graph
  // vertex properties
  property<vertex_name_t, string,
     property <vertex_index_t, int>  //property added by Nicolas Wicker
     >,
  // edge properties
  property<edge_weight_t, double,
     property <edge_name_t, string>
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
  for(tie(vitr,vend)=vertices(g); vitr != vend; ++vitr){
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

namespace std {
  template <typename T>
  istream& operator >> (istream& in, pair<T,T>& p)
  {
    in >> p.first >> p.second;
    return in;
  }
}

Graph& get_graph(pair<Graph,VertexPair>& _data){
  return _data.first;
}

Vertex& get_start_vertex(pair<Graph,VertexPair>& _data){
  return _data.second.first;
}

Vertex& get_end_vertex(pair<Graph,VertexPair>& _data){
  return _data.second.second;
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

pair<Graph,VertexPair> SCOPE_R_get_st_graph_from(SEXP nodes, SEXP edges,SEXP edge_weights,int *nt);

/*********************************************************************************/
/*                                                                               */
/*   Nicolas Wicker addition                                                     */
/*                                                                               */
/*********************************************************************************/

/************************************************************************/
/*                                                                      */
/*Procedure compute empirical score distribution for a given path length*/
/*by random sampling, the sampling is done by Metropolis algorithm      */
/*                                                                      */
/************************************************************************/

void computeRandomScores(int maximumPathLength,const Vertex s,Graph &g,int nbSamples,
                         int nbWarmingSteps,double **randomScores)
{
  //int currentIndex,currentLength,precedingVertexIndex,deadEnd,pvalue=0;
  int i,j,k,length,noLoop,vindex;
  int iterator,randomValue,nbValidNeighbours;
  EdgeWeightMap  weight = get(edge_weight, g);
  AdjacencyIter aitr, aend;
  Vertex currentVertex;
  VertexIndexMap vertexIndexMap=get(vertex_index, g);
  int *indices,validVertexFound,validPath;
  int nbVertices;

  int warmingTime;
  double logCurrentPathProbability,logProposalPathProbability;
  double currentScore,proposalScore;

  indices=new int[maximumPathLength+1];
  /*vertices array to obtain a vertex from an index*/
  nbVertices = num_vertices(g);
  Vertex *verticesArray=new Vertex[nbVertices];
  VertexIter vitr,vend;
  for(tie(vitr,vend)=vertices(g); vitr != vend; ++vitr) {
    verticesArray[vertexIndexMap[*vitr]]=*vitr;
  }

  for(length=1;length<=maximumPathLength;length++) {
    int nbChanges=0;
    currentScore=0;
    logCurrentPathProbability=DBL_MAX;
    int nbFailures=0;
    int currentSampleSize=0;
    for(iterator=1;iterator<=nbWarmingSteps*nbSamples;iterator++) {
        if((iterator%nbWarmingSteps)==0) {
            warmingTime=0;
        } else {
            warmingTime=1;
        }
        validPath=0;
        while(validPath==0) {
            /*here the paths building procedure starts*/
            logProposalPathProbability=0.0;
            validPath=1;
            proposalScore=0.0;
            currentVertex=verticesArray[(int)floor(runif(0,nbVertices))];
            indices[0]=vertexIndexMap[currentVertex];
            for(i=0;i<length;i++) {
                if(out_degree(currentVertex,g)==0) {
                    break;
                }
                /*check how many non-visited neighbours there are*/
                nbValidNeighbours=0;
                for(tie(aitr,aend)=adjacent_vertices(currentVertex,g); aitr != aend; ++aitr) {
                    /*must verify first if there is no loop*/
                    noLoop=1;
                    vindex=vertexIndexMap[*aitr];
                    for(k=i;k>=0;k--) {
                        if(vindex==indices[k]) {
                            noLoop=0;
                            break;
                        }
                    }
                    if(noLoop==1) {
                        nbValidNeighbours++;
                    }
                }
                if(nbValidNeighbours>0) {
                    logProposalPathProbability-=log((double)nbValidNeighbours);
                    validVertexFound=0;
                    while(validVertexFound==0) {
                        randomValue=(int)floor(runif(0,out_degree(currentVertex,g)));
                        j=0;
                        for(tie(aitr,aend)=adjacent_vertices(currentVertex,g); aitr != aend; ++aitr) {
                            if(j==randomValue) {
                                proposalScore+=(double)weight[edge(currentVertex,*aitr,g).first];
                                /*must verify first if there is no loop*/
                                noLoop=1;
                                vindex=vertexIndexMap[*aitr];
                                for(k=i;k>=0;k--) {
                                    if(vindex==indices[k]) {
                                        noLoop=0;
                                        break;
                                    }
                                }
                                if(noLoop==1) {
                                    currentVertex=*aitr;
                                    indices[i+1]=vindex;
                                    validVertexFound=1;
                                }
                                break;
                            }
                            j++;
                        }
                    }
                } else {
                    /*the random path has path reached a dead-end*/
                    validPath=0;
                    nbFailures++;
                    break;
                }
            }
        }
        /*end of the path building*/
        if(nbFailures<currentSampleSize) {
            logProposalPathProbability-=log(1.0-(double)nbFailures/(double)currentSampleSize);
        }
        if(runif(0,1)<exp(logCurrentPathProbability-logProposalPathProbability)) {
            nbChanges++;
            currentScore=proposalScore;
            logCurrentPathProbability=logProposalPathProbability;
        }
        if(warmingTime==0) {
            randomScores[length][currentSampleSize]=currentScore;
            logCurrentPathProbability=DBL_MAX;
            currentSampleSize++;
        }
     }
  }
  delete []indices;
  delete []verticesArray;
}

/************************************************************************/
/*                                                                      */
/*Procedure compute empirical score distribution for a given path length*/
/*by random sampling, NOT Metropolis sampling							*/
/*Ahmed.																*/
/*                                                                      */
/************************************************************************/

void computeRandomScoresRandomSampling(int maxlength,int nbEdges, int nbIterations,double *weights, double** randomScores) {

  int i,j;
  int iterator,randomValue;

  //double randomScore;
  Rprintf("sampling started");
  for(i=0;i<maxlength;i++) {
	  randomScores[i] = new double[nbIterations];
	for(iterator=0;iterator<nbIterations;iterator++) {
		randomScores[i][iterator] = 0.0;
    	for(j=0; j<i+1; j++){
    		randomValue=(int)floor(runif(0,nbEdges-1));
    		randomScores[i][iterator]+=weights[randomValue];
    	}
    }
  }
  Rprintf("sampling finished");
}

/******************************************************/
/*                                                    */
/*Procedure to compute pvalues by path random sampling*/
/*the sampling is done by Metropolis algorithm        */
/*                                                    */
/******************************************************/

double computePvalue(double score,int length,int nbSamples,double **randomScores) {
  int a,b,c;
  //double pvalue;

  a=0;
  b=nbSamples-1;
  if(randomScores[length][0]>=score) {
    return 0;
  }
  while(true) {
    c=(a+b)/2;
    if(randomScores[length][c]>=score) {
      b=c;
    } else {
      a=c;
    }
    if(b==(a+1)) {
      break;
    }
  }
  return a/(double)nbSamples;
}

/***************************************************************************/
/*                                                                         */
/*Procedure to compute pvalues by weights random sampling (does not work !)*/
/*                                                                         */
/***************************************************************************/

double computePvalue2(double score,int length,int nbEdges,double *weights) {

  int i;
  int iterator,nbIterations,randomValue;

  nbIterations=100000;
  int pvalue=0;
  double randomScore;
  for(iterator=0;iterator<nbIterations;iterator++) {
    randomScore=0.0;
    for(i=0;i<length;i++) {
      randomValue=(int)floor(runif(0,nbEdges-1));
      randomScore+=weights[randomValue];
    }
    if(randomScore<=score) {
      pvalue++;
    }
  }
  return (double)pvalue/(double)nbIterations;
}

/************************************************/
/*                                              */
/*Procedure to find out the shortest pvalue path*/
/*                                              */
/************************************************/

//struct st_path_with_deviation
SEXP st_shortestPvalue_path(const Vertex s,const Vertex t,Graph &g,
                       int maximumPathLength,int nbSamples,double **randomScores,double alpha) {
  int i,j,n;
  int sindex,tindex;
  int length,uindex,vindex,noLoop,precedingVertexIndex,currentIndex,currentLength;
  Vertex u,v;
  double rhs;
  int **pathsMatrix; /*pathsMatrix indicates for each pari (vertexIndex,length) the preceding vertexIndex
                           for the shortest path of length l*/
  double **scoresMatrix;
  AdjacencyIter aitr, aend;
  n = num_vertices(g);

  /*call shortest path function*/
  struct st_path_with_deviation shortestPath=st_shortest_path(s,t,g);

  VertexIndexMap vertexIndexMap=get(vertex_index, g);
  EdgeWeightMap  weight_map = get(edge_weight, g);
  VertexNameMap vertexNameMap = get(vertex_name, g);
  EdgeNameMap   edge_name_map   = get(edge_name, g);

  /*memory allocation*/
  /*vertices array to obtain a vertex from an index*/
  Vertex *verticesArray=new Vertex[n];
  scoresMatrix=new double *[n];
  pathsMatrix=new int *[n];
  /*end of memory allocation*/

   /*initialization*/
  for(i=0;i<n;i++) {
    scoresMatrix[i]=new double[n];
    pathsMatrix[i]=new int[n];
  }
  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      scoresMatrix[i][j]=-1.0;
    }
  }
  VertexIter vitr,vend;
  for(tie(vitr,vend)=vertices(g); vitr != vend; ++vitr) {
    verticesArray[vertexIndexMap[*vitr]]=*vitr;
  }

  sindex=vertexIndexMap[s];
  tindex=vertexIndexMap[t];
  scoresMatrix[sindex][0]=0.0;

  for(length=0;length<maximumPathLength;length++) {
    for(tie(vitr,vend)=vertices(g); vitr != vend; ++vitr) {
        u=*vitr;
        uindex=vertexIndexMap[u];

        /*impossible value*/

        if(scoresMatrix[uindex][length]>-0.5) {

            for(tie(aitr,aend)=adjacent_vertices(u,g); aitr != aend; ++aitr) {
                v   = *aitr;
                vindex=vertexIndexMap[v];
                rhs = scoresMatrix[uindex][length]+(double)weight_map[edge(u,v,g).first];

                /*must verify first if there is no loop*/
                noLoop=1;
                currentIndex=uindex;
                currentLength=length;

                while(currentIndex!=sindex) {

                    precedingVertexIndex=pathsMatrix[currentIndex][currentLength];

                    if(precedingVertexIndex==vindex) {
                        noLoop=0;
                        break;
                    } else {
                        currentIndex=precedingVertexIndex;
                        currentLength--;
                    }
                }

                if ((noLoop==1)&&((scoresMatrix[vindex][length+1] > rhs) || (scoresMatrix[vindex][length+1]<-0.5)))
                {
                    scoresMatrix[vindex][length+1] = rhs;
                    pathsMatrix[vindex][length+1] = uindex;

                }
           }
        }
     }
  }

  /*generate pvalues*/
  double *pvalues;
  EdgeIter eitr,eend;
  double *weights;
  SizeType nbEdges;

  nbEdges=num_edges(g);
  weights=new double[nbEdges];

  i=0;
  for(tie(eitr,eend)=edges(g); eitr != eend; ++eitr) {
    weights[i]=weight_map[*eitr];
    i++;
  }

  pvalues=new double[n];
  /*i below is the number of edges, so that the path length is equal to i+1*/
  /*the first edge is sorted in pathsMatrix[vindex][1] */
  deque <Vertex>path;
  deque <double>pathWeights;
  Edge e_temp;

  SEXP genes = R_NilValue,compounds = R_NilValue,pweights = R_NilValue,distance = R_NilValue,pp = R_NilValue;

  int allocated = 0;

  for(length=0;length<maximumPathLength;length++) {
	//cout<<length << ":"<<scoresMatrix[tindex][length] << "  ";
    if(scoresMatrix[tindex][length]>0) {
            pvalues[length]=computePvalue(scoresMatrix[tindex][length],length,nbSamples,randomScores);
            //cout << "scores: "<< scoresMatrix[tindex][length]<< ",pval: " << pvalues[length] <<endl;

            if(pvalues[length] < alpha) {

                Rf_protect(compounds = Rf_allocVector(STRSXP,length-1));
                Rf_protect(genes = Rf_allocVector(STRSXP,length));
                Rf_protect(pweights = Rf_allocVector(REALSXP,length-1));
                Rf_protect(distance = Rf_allocVector(REALSXP,1));
                Rf_protect(pp =  Rf_allocVector(REALSXP,1));

                REAL(distance)[0] = scoresMatrix[tindex][length];
                REAL(pp)[0] = pvalues[length];
                allocated = 1;

                /*path display*/
                currentIndex=tindex;
                currentLength=length;
                path.clear();
                while(currentIndex!=sindex) {
                    precedingVertexIndex=pathsMatrix[currentIndex][currentLength];
                    u=verticesArray[precedingVertexIndex];
                    path.push_front(u);
                    v=verticesArray[currentIndex];
                    currentIndex=precedingVertexIndex;
                    currentLength--;
                }
                //cout << "length: " <<length<<"size: " << path.size() <<endl;


                //u=path.front();
                //SET_STRING_ELT(compounds,0,Rf_mkChar(vertexNameMap[u].c_str()));

                //v=u;
                path.pop_front();
                pathWeights.clear();

                u=path.front();
                SET_STRING_ELT(genes,0,Rf_mkChar(vertexNameMap[u].c_str()));
                path.pop_front();
                v=u;

                int i = 0;
                while(not path.empty()) {
                    u=path.front();
                    e_temp = edge(v,u,g).first;
                    pathWeights.push_back(weight_map[e_temp]);
                    v=u;

                    SET_STRING_ELT(compounds,i,Rf_mkChar(edge_name_map[e_temp].c_str()));
                    SET_STRING_ELT(genes,i+1,Rf_mkChar(vertexNameMap[u].c_str()));

                    i = i + 1;

                    path.pop_front();
                }
                e_temp = edge(v,t,g).first;
                pathWeights.push_back(weight_map[e_temp]);

                SET_STRING_ELT(compounds,length-2,Rf_mkChar(edge_name_map[e_temp].c_str()));
                SET_STRING_ELT(genes,length-1,Rf_mkChar(vertexNameMap[t].c_str()));

                i = 0;
                while(not pathWeights.empty()) {
                    REAL(pweights)[i] = pathWeights.front();
                    pathWeights.pop_front();
                    i = i + 1;
                }

                break;
                /*end of path display*/

        } else if(pvalues[length]>0.1) {
               break;
        }
    }
  }

  /*memory desallocation*/
  for(i=0;i<n;i++) {
    delete [] scoresMatrix[i];
    delete [] pathsMatrix[i];
  }
  delete [] scoresMatrix;
  delete [] pathsMatrix;
  delete [] pvalues;
  delete [] weights;
  /*end of memory desallocation*/

  SEXP np, dimnames;
  if (allocated == 1) {
    Rf_protect(np = Rf_allocVector(VECSXP, 5));
    Rf_protect(dimnames = Rf_allocVector(VECSXP,5));
    SET_VECTOR_ELT(np,0,genes);
    SET_VECTOR_ELT(dimnames,0,Rf_mkString("genes"));
    SET_VECTOR_ELT(np,1,compounds);
    SET_VECTOR_ELT(dimnames,1,Rf_mkString("compounds"));
    SET_VECTOR_ELT(np,2,pweights);
    SET_VECTOR_ELT(dimnames,2,Rf_mkString("weights"));
    SET_VECTOR_ELT(np,3,distance);
    SET_VECTOR_ELT(dimnames,3,Rf_mkString("distance"));
    SET_VECTOR_ELT(np,4,pp);
    SET_VECTOR_ELT(dimnames,4,Rf_mkString("pvalue"));
    Rf_setAttrib(np,R_NamesSymbol,dimnames);
    Rf_unprotect(7);
    return(np);
  } else {
    return R_NilValue;
  }

}

/*********************************************************************************/
/*                                                                               */
/*   End of Nicolas Wicker addition                                              */
/*                                                                               */
/*********************************************************************************/

extern "C" SEXP scope(SEXP node_list,
                      SEXP edge_list,
                      SEXP edge_weights,
                      SEXP SAMPLEDPATHS,
                      SEXP ALPHA,
                      SEXP ECHO)
{
  Vertex s,t;
  Graph g;
  int i,start,end,counter;
  int echo = INTEGER(ECHO)[0];
  double **randomScores;
  int maxSamplePathLength = Rf_ncols(SAMPLEDPATHS), numberPathSamples = Rf_nrows(SAMPLEDPATHS);
  double alpha = REAL(ALPHA)[0];

  int nt = 0;
  pair<Graph,VertexPair> data = SCOPE_R_get_st_graph_from(node_list, edge_list, edge_weights,&nt);


  g = get_graph(data);
  s = get_start_vertex(data);
  t = get_end_vertex(data);

  // if sampled paths were not provided, use random edge sampling.
  if(SAMPLEDPATHS!=R_NilValue){
	  double *rs = REAL(SAMPLEDPATHS);
	  randomScores= new double *[maxSamplePathLength];
	  for (i = 0;i < maxSamplePathLength;i = i + 1) randomScores[i] = &rs[i*numberPathSamples];
  }else{
	  maxSamplePathLength = num_vertices(g);
	  numberPathSamples = 10000;
	  randomScores= new double *[maxSamplePathLength];
	  computeRandomScoresRandomSampling(maxSamplePathLength, num_edges(g), numberPathSamples,REAL(edge_weights), randomScores);
  }

  if(s == numeric_limits<Vertex>::max() or t == numeric_limits<Vertex>::max()){
    Rprintf("No vertex start or end vertices found.");
    return R_NilValue;
  }

  start=0;
  end=num_vertices(g);

  VertexNameMap vertex_name_map = get(vertex_name, g);
  VertexIter vitr,vend;
  AdjacencyIter aitr, aend;

  //struct st_path_with_deviation p, q;
  SEXP newpath,allpaths;
  Rf_protect(allpaths = Rf_allocVector(VECSXP,nt));
  vector<string> reachedTargets;
  vector<string>::iterator it;

  counter=0;
  if (echo == 1) Rprintf("There are %d nodes in the neighborhood", nt);
  int vectid = 0;
  for(tie(vitr,vend)=vertices(g); vitr != vend; ++vitr) {
    for(tie(aitr,aend)=adjacent_vertices(*vitr,g); aitr != aend; ++aitr) {
        if(*aitr==t) {
            string target(vertex_name_map[*vitr]);
            string last_compound = target;

            if (echo == 1) Rprintf("TARGET : %s %d/ %d", last_compound.c_str() ,vectid+1 ,nt);
            int alreadyreached = 0;
            for (it=reachedTargets.begin();it < reachedTargets.end();it++) {
                if (last_compound.compare(*it) == 0) {
                    alreadyreached = 1;
                    if (echo == 1) Rprintf(" - Already found a path to %s\n", last_compound.c_str());
                }
            }
            if (alreadyreached == 0) {
                if((start<=counter)&&(counter<end)) {
                    newpath = st_shortestPvalue_path(s,*vitr,g,maxSamplePathLength,numberPathSamples,randomScores,alpha);
                    if (newpath != R_NilValue) {
                        SET_VECTOR_ELT(allpaths,vectid,newpath);
                        reachedTargets.push_back(last_compound);
                         if (echo == 1) Rprintf(" - Found a path to %s\n", last_compound.c_str());;
                    } else {
                         if (echo == 1) Rprintf(" - Node out of scope\n");
                    }
                }
            }
            vectid = vectid + 1;
        }
    }
    counter++;
  }

  SEXP scopecompounds, out, dimnames;

  Rf_protect(scopecompounds = Rf_allocVector(STRSXP,reachedTargets.size()));
  for (int s = 0;s < (int)reachedTargets.size();s++) {
      SET_STRING_ELT(scopecompounds,s,Rf_mkChar(reachedTargets.at(s).c_str()));
  }

  Rf_protect(out = Rf_allocVector(VECSXP,2));
  Rf_protect(dimnames = Rf_allocVector(VECSXP,2));
  SET_VECTOR_ELT(out,0,allpaths);
  SET_VECTOR_ELT(dimnames,0,Rf_mkString("paths"));
  SET_VECTOR_ELT(out,1,scopecompounds);
  SET_VECTOR_ELT(dimnames,1,Rf_mkString("scope"));
  Rf_setAttrib(out,R_NamesSymbol,dimnames);

  /*memory desallocation*/
  delete [] randomScores;
  Rf_unprotect(4);

  return(out);
}

pair<Graph,VertexPair> SCOPE_R_get_st_graph_from(SEXP nodes, SEXP edges,SEXP edge_weights,int *nt) {

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
  string      edge_label;

  EdgeNameMap   edge_name_map   = get(edge_name, g);

  // define R types
  from_idxs = VECTOR_ELT(edges,0);
  to_idxs = VECTOR_ELT(edges,1);
  labels = VECTOR_ELT(edges,2);

  *nt = 0;
  for (int i = 0;i < LENGTH(from_idxs);i = i + 1) {
	v_pair.first = (int)(REAL(from_idxs)[i] - 1);
	v_pair.second = (int)(REAL(to_idxs)[i] - 1);
	if ((int)(REAL(to_idxs)[i] - 1) == (int)t) *nt = *nt + 1;

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

extern "C" SEXP samplepaths(SEXP node_list,
                      SEXP edge_list,
                      SEXP edge_weights,
                      SEXP MAXPATHLENGTH,
                      SEXP SAMPLEPATHS,
                      SEXP WARMUPSTEPS)
{
  Vertex s;
  Graph g;

  int maximumPathLength=INTEGER(MAXPATHLENGTH)[0];
  int nbSamples=INTEGER(SAMPLEPATHS)[0];
  int nbWarmingSteps=INTEGER(WARMUPSTEPS)[0];
  int i,j;

  SEXP SCORES;
  Rf_protect(SCORES = Rf_allocVector(REALSXP,nbSamples*(maximumPathLength+1)));
  double *rs = REAL(SCORES);

  double **randomScores= new double *[(maximumPathLength+1)];
  for (i = 0;i < (maximumPathLength+1);i = i + 1) randomScores[i] = &rs[i*(nbSamples)];

  int nt = 0;
  pair<Graph,VertexPair> data = SCOPE_R_get_st_graph_from(node_list, edge_list, edge_weights,&nt);

  g = get_graph(data);
  s = get_start_vertex(data);

  /*compute random scores for length ranging from 1 to maximumPathLength*/
  GetRNGstate();
  computeRandomScores(maximumPathLength,s,g,nbSamples,nbWarmingSteps,randomScores);
  PutRNGstate();
  std::vector<double> orderedScores;
  for(i=1;i<=maximumPathLength;i++) {
    orderedScores.clear();
    for(j=0;j<nbSamples;j++) {
        orderedScores.push_back(randomScores[i][j]);
    }
    std::sort(orderedScores.begin(),orderedScores.end());
    for(j=0;j<nbSamples;j++) {
        randomScores[i][j]=orderedScores[j];
    }
  }

  Rf_unprotect(1);
  delete [] randomScores;

  return(SCORES);
}
