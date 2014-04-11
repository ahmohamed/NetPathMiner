#include "init.h"


template <class T>
typename vector<T>::size_type elem_pos(vector<T> v, const T &e){
	return( find(v.begin(), v.end(), e) - v.begin() );
}

template <class T>
typename vector<T>::size_type add_elem(vector<T> &v,const T &e){
	typename vector<T>::size_type pos = elem_pos(v,e);
	if(pos == v.size())
		v.push_back(e);

	return pos;
}

int compare (const void * a, const void * b)
{
  return ( *(double*)a - *(double*)b );
}
double median(double x[], int n) {
    if(n==0){return(NA_REAL);}
    if(n==1){return(x[0]);}


    qsort(x, n, sizeof(x[0]), compare);
    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((x[n/2] + x[n/2 - 1]) / 2.0);
    } else {
        // else return the element in the middle
        return x[n/2];
    }
}

void corEdgeWeights(double * X,
    int * EDGELIST,
    int * SAMEGENE,
    double * WEIGHT,
    int *NEDGES,
    int * NOBS,
    int * NCOR)
{
    int nobs = (int)(*NOBS);
    int nedges = (int)(*NEDGES);
    const int ncor = (int)(*NCOR);
    int indx, i;

    // For each edge
    for (indx = 0;indx < nedges;indx = indx + 1) {
        int to_indx = EDGELIST[indx+nedges];
        int from_indx = EDGELIST[indx];

        if(to_indx == NA_INTEGER || from_indx == NA_INTEGER){
        	WEIGHT[indx] = NA_REAL;
        	continue;
        }

        WEIGHT[indx] = 0.0; // if all else fails the weight will be assigned to -1 i.e. the most inprobable edge

        // Compute the correlation
        if (SAMEGENE[indx] == 0) {
        	double* corlist = new double[ncor];

            for(int j=0; j<ncor; j++){
				double Exy = 0.0, Exx = 0.0, Ex = 0.0, Eyy = 0.0, Ey = 0.0;
				double xp = 0.0, yp = 0.0;
				double n = (double)nobs;

				for (i = 0;i < nobs;i = i + 1) {
					if(ncor >1){  //If multiple correlations, sample the columns and take the median.
						int sample = rand() % nobs;
						xp = X[from_indx*nobs + sample]; yp = X[to_indx*nobs + sample];
					}else{
					xp = X[from_indx*nobs + i]; yp = X[to_indx*nobs + i];
					}

					if (!isnan(xp) && !isnan(yp)) {
						Ex = Ex + xp; Exx = Exx + xp*xp; Ey = Ey + yp; Eyy = Eyy + yp*yp; Exy = Exy + xp*yp;
					} else n = n - 1.0; // If it is a missing value skip that observation completely and reduce the dataset size by 1.
				}
				if (n > 2) {
					if (Exy != 0.0 && Exx != 0.0 && Eyy != 0.0 && Ex != 0.0 && Ey != 0.0)
						corlist[j] = (n*Exy - Ex*Ey)/ sqrt( (n*Exx - Ex*Ex) * (n*Eyy - Ey*Ey) );
				}
            }
            WEIGHT[indx] = median(corlist, ncor);
            delete[] corlist;
        } else {
            // If it is the same gene set to minimum of -1.0 (penalty)
            WEIGHT[indx] = -1.0;
        }
    }
}

SEXP expand_complexes(SEXP ATTR_LS, SEXP EL, SEXP V, SEXP EXPAND, SEXP MISSING){
	/* Processing arguments */
	vector<vector<string> > attr_ls;
	vector<string> v_name;
	bool duplicate = ((string)CHAR(STRING_ELT(EXPAND,0)))=="duplicate";
	string missing = CHAR(STRING_ELT(MISSING,0));

	for(int i=0; i<LENGTH(ATTR_LS); i++){
		SEXP VEC = AS_CHARACTER(VECTOR_ELT(ATTR_LS, i));
		attr_ls.push_back(vector<string>());
		if(LENGTH(VEC)==0)
			continue;
		for(int j=0;j<LENGTH(VEC);j++)
			attr_ls[i].push_back(CHAR(STRING_ELT(VEC,j)));
	}

	for(int i=0;i<LENGTH(V);i++)
		v_name.push_back( CHAR(STRING_ELT(V,i)) );

	/* Holder vectors */
	vector<string> vertices; //Vertex list
	vector<int> edges;		//Edge list
	vector< vector<int> > parents;	//For each expanded vertex, keep its parent(s) indices, for attribute inheritance.
	vector<int> edge_parents;
	vector<int> non_gene; //Keep track of unexpandable vertices (missing annotation), to remove them later.

	for(int i=0; i<LENGTH(EL); i+=2){
		vector<size_t> el1_pos, el2_pos;
		int el1,el2;
		el1 = INTEGER(EL)[i]; el2 = INTEGER(EL)[i+1];

		//Add annotations as vertices, store their indices in el_pos vectors
		for(size_t j= 0; j< attr_ls[el1].size(); j++){
			string el_name;
			if(!duplicate)
				el_name = attr_ls[el1][j];
			else{el_name = attr_ls[el1][j] + ((string)"##") + v_name[el1];}

			el1_pos.push_back(add_elem(vertices, el_name));
			if(el1_pos[j]==parents.size())
				parents.push_back(vector<int>());
			add_elem(parents[el1_pos[j]],el1);
		}

		for(size_t j=0; j<attr_ls[el2].size(); j++){
			string el_name;
			if(!duplicate)
				el_name = attr_ls[el2][j];
			else{el_name = attr_ls[el2][j] + ((string)"##") + v_name[el2];}

			el2_pos.push_back(add_elem(vertices, el_name));
			if(el2_pos[j]==parents.size())
				parents.push_back(vector<int>());
			add_elem(parents[el2_pos[j]],el2);
		}

		/* *If missing attribute vertices are to be removed,no edges will be added.
		 * *If they are kept, the vertex original ID (rather than its annotation) is added
		 * 	as vertex in the new graph.
		 * *If they are to be removed, and their neighbours connected, indices of them
		 *  are kept, for later use (in R counterpart).
		 */

		for(vector<int>::size_type j = 0; j != el1_pos.size(); j++){
			for(vector<int>::size_type k = 0; k != el2_pos.size(); k++){
				edges.push_back(el1_pos[j]);
				edges.push_back(el2_pos[k]);
				edge_parents.push_back( (int)(i/2 + 1) );
			}
		}
	}//looping over the edgelist

	/* Storing the results in R objects */
	SEXP OUT, NAMES, VERTICES, EDGES, E_PARENTS, PARENTS, RECONNECT;
	PROTECT( VERTICES = NEW_STRING(vertices.size()) );
	PROTECT( EDGES = NEW_INTEGER(edges.size()) );
	PROTECT( E_PARENTS = NEW_INTEGER(edge_parents.size()) );
	PROTECT( RECONNECT = NEW_INTEGER(non_gene.size()) );
	PROTECT( PARENTS = NEW_LIST(parents.size()) );

	for(size_t i=0; i<vertices.size(); i++)
		SET_STRING_ELT(VERTICES, i, mkChar( vertices[i].c_str() ));

	for(size_t i=0; i<edges.size(); i++)
		INTEGER(EDGES)[i] = edges[i]+1;

	for(size_t i=0; i<edge_parents.size(); i++)
		INTEGER(E_PARENTS)[i] = edge_parents[i];

	for(size_t i=0; i<non_gene.size(); i++)
		INTEGER(RECONNECT)[i] = non_gene[i]+1;

	for(size_t i=0; i<parents.size(); i++){
		SEXP PARENT_i;
		PROTECT(PARENT_i = NEW_INTEGER(parents[i].size()));
		for(size_t j=0; j<parents[i].size(); j++)
			INTEGER(PARENT_i)[j] = parents[i][j]+1;

		SET_VECTOR_ELT(PARENTS, i, PARENT_i);
		UNPROTECT(1);
	}

	PROTECT( OUT = NEW_LIST(5));
	PROTECT( NAMES = NEW_STRING(5));
	SET_VECTOR_ELT(OUT, 0, VERTICES); SET_STRING_ELT(NAMES, 0, mkChar("vertices"));
	SET_VECTOR_ELT(OUT, 1, EDGES);	SET_STRING_ELT(NAMES, 1, mkChar("edges"));
	SET_VECTOR_ELT(OUT, 2, RECONNECT);	SET_STRING_ELT(NAMES, 2, mkChar("reconnect"));
	SET_VECTOR_ELT(OUT, 3, PARENTS);	SET_STRING_ELT(NAMES, 3, mkChar("parents"));
	SET_VECTOR_ELT(OUT, 4, E_PARENTS);	SET_STRING_ELT(NAMES, 4, mkChar("e.parents"));

	setAttrib(OUT,R_NamesSymbol,NAMES);
	UNPROTECT(7);
	return(OUT);
}
