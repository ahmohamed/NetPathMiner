#ifdef HAVE_SBML
#ifndef __sbml_interface__h_
#define __sbml_interface__h_

#include <sstream>
#include <sbml/SBMLTypes.h>
#include "handlesegfault.h"
#include "init.h"

SEXP getReactionList(Model *model, const vector<string> &attr_terms, vector<string> &species, bool verbose);
SEXP getSpeciesFrame(Model *model, vector<string> species, const vector<string> &attr_terms);

SEXP get_species_info(Model *model, const string species, const vector<string> &attr_terms);
void readsbml_sign_int(Model *model, vector<string> &species, vector<size_t> &non_gene,
						vector<SEXP> &info, vector<size_t> &edges,
						const vector<string> &attr_terms, bool verbose);

void get_MIRIAM(XMLNode* rdf, const vector<string> &terms, vector< vector<string> > &values, vector<string> &names);
template <class T> size_t elem_pos(vector<T> v, const T &e);
template <class T> size_t add_elem(vector<T> &v, const T &e);
const char* URL_decode(char* URL);
bool not_alnum(char c);

#undef length
#endif
#endif
