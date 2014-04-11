#ifdef HAVE_XML
#include <libxml/xmlreader.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <sstream>

#include "handlesegfault.h"
#include "init.h"

/* Declaration of functions */
void readkgml_sign_int(const char* filename, vector<string> &vertices,
						vector<int> &edges,	vector< vector<string> > &attr,
						vector< vector<string> > &pathway_attr, bool expand_complexes,
						bool verbose);
xmlNodePtr node_by_id(char* id, const char* type, xmlXPathContextPtr &xpathCtx);
xmlNodePtr node_by_attr_val(const char* attr, char* val, const char* type, xmlXPathContextPtr &xpathCtx);
char* get_attr(xmlNodePtr node,const char* attr_name);
char* attr_by_id(char* id, const char* attr_name, xmlXPathContextPtr &xpathCtx);
std::vector<std::string> split(const std::string &s, char delim);
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
char* get_group_components(char* id, xmlXPathContextPtr &xpathCtx);

template <class T>
size_t elem_pos(vector<T> v, T &e);

template <class T>
bool elem_in_vector(vector<T> v, T &e);



SEXP readkgmlfile(SEXP FILENAME, SEXP VERBOSE) {
	handle_segfault_KGML();

	const char *filename = CHAR(STRING_ELT(FILENAME,0));
	bool verbose = LOGICAL(VERBOSE)[0];

	xmlDocPtr doc;
	xmlXPathContextPtr xpathCtx;
	xmlXPathObjectPtr nodes;


	if(verbose)	Rprintf("Processing KGML file: %s",filename);

	/* Load XML document */
	doc = xmlParseFile(filename);
	if (doc == NULL) {
		Rf_warningcall(mkChar(filename), "Unable to parse file");
		if(verbose)	Rprintf(": Error.\n");
		return(R_NilValue);
	}

	/* Check it is a kegg pathway file */
	if(doc->intSubset == NULL ||
	   strcmp( (char *) (doc->intSubset->name), "pathway") != 0 )
	   //strncmp( (char *) (doc->intSubset->SystemID), "http://www.kegg.jp/kegg/", 24) !=0)
	{
		Rf_warningcall(mkChar(filename), "File is not KEGG pathway file");
		xmlFreeDoc(doc);
		if(verbose)	Rprintf(": Error.\n");
		return(R_NilValue);
	}

	/* Get pathway information :*/
	xmlNodePtr pathway =  xmlDocGetRootElement(doc);
	if(pathway == NULL || strcmp( (char *) (pathway->name), "pathway") != 0){
		Rf_warningcall(mkChar(filename), "No pathways in file");
		xmlFreeDoc(doc);
		if(verbose)	Rprintf(": Error.\n");
		return(R_NilValue);
	}
	const char* pathwayId = get_attr(pathway, "name");
	if(!pathwayId){
		Rf_warningcall(mkChar(filename), "Pathway ID not found in file. Using file name instead.");
		pathwayId = filename;
	}else{
		pathwayId +=5; //Remove "path:" leading characters//
	}

	const char* pathwayTitle = get_attr(pathway, "title");
	if(!pathwayTitle){
		Rf_warningcall(mkChar(pathwayId), "Pathway title not found in file.");
		pathwayTitle = "";
	}
	if(verbose)	Rprintf(" \"%s\"",pathwayTitle);

	/* Create xpath evaluation context */
	xpathCtx = xmlXPathNewContext(doc);
	if(xpathCtx == NULL) {
		Rf_warningcall(mkChar(filename), "Unable to create new XPath context");
		xmlFreeDoc(doc);
		if(verbose)	Rprintf(": Error.\n");
		return(R_NilValue);
	}

	/* Evaluate xpath expression */
	nodes = xmlXPathEvalExpression((xmlChar *) "//reaction", xpathCtx);
	if(nodes == NULL || nodes->nodesetval == NULL || nodes->nodesetval->nodeNr == 0) {
		Rf_warningcall(mkChar(pathwayId), "Pathway contains no reactions");
		xmlXPathFreeContext(xpathCtx);
		xmlFreeDoc(doc);
		if(verbose)	Rprintf(": Error.\n");
		return(R_NilValue);
	}


	/* Parse XML Reactions*/
	xmlNodePtr curReaction;
    int size;
    int i;
    size = (nodes->nodesetval) ? nodes->nodesetval->nodeNr : 0;

    SEXP REACTIONLIST,ID;
    PROTECT(REACTIONLIST = allocVector(VECSXP,size));
    PROTECT(ID = allocVector(STRSXP,size));

    if(verbose)	Rprintf(": %d reactions found.\n",size);

    char* temp; //template for string processing
    for(i = 0; i < size; ++i) {
    	curReaction = nodes->nodesetval->nodeTab[i];
    	xpathCtx->node = curReaction;

    	const char *name = get_attr(curReaction, "name");
    	SET_STRING_ELT(ID,i,mkChar(name));

    	SEXP REACTION,REACTIONNAMES;
		SEXP NAME, REVERSIBLE, REACTANTS, RSTOIC, PRODUCTS, PSTOIC, GENES, PATHWAY,
		KEGG_PATHWAY, KEGG_REACTION, KEGG_GENES, NCBI_GENE;

		PROTECT(REACTION = allocVector(VECSXP,13));
		PROTECT(REACTIONNAMES = allocVector(STRSXP,13));


		PROTECT(NAME = allocVector(STRSXP,1));
		SET_STRING_ELT(NAME,0, mkChar(name));
		SET_VECTOR_ELT(REACTION,0,NAME); SET_STRING_ELT(REACTIONNAMES,0,mkChar("name"));

		PROTECT(REVERSIBLE = allocVector(LGLSXP,1));
		LOGICAL(REVERSIBLE)[0] = strcmp(get_attr(curReaction, "type"), "irreversible") != 0;
		SET_VECTOR_ELT(REACTION,1,REVERSIBLE); SET_STRING_ELT(REACTIONNAMES,1,mkChar("reversible"));
		//cout << (char *) xmlGetProp(curReaction,(const xmlChar *)"type") <<endl;


		xmlNodeSetPtr reactantNodes = xmlXPathEvalExpression( (const xmlChar *) "./substrate", xpathCtx)->nodesetval;
		int numOfReactants = (reactantNodes) ? reactantNodes->nodeNr : 0;
		PROTECT(REACTANTS = allocVector(STRSXP,numOfReactants));
		PROTECT(RSTOIC = allocVector(REALSXP,numOfReactants));

		for (int r = 0;r < numOfReactants;r++) {
			temp = get_attr(reactantNodes->nodeTab[r], "name");
			SET_STRING_ELT(REACTANTS,r,mkChar(temp+4));
			REAL(RSTOIC)[r] = NA_REAL;
		}
		SET_VECTOR_ELT(REACTION,2,REACTANTS); SET_STRING_ELT(REACTIONNAMES,2,mkChar("reactants"));
		SET_VECTOR_ELT(REACTION,3,RSTOIC);SET_STRING_ELT(REACTIONNAMES,3,mkChar("reactant.stoichiometry"));
		//cout << "Reactant level :|" << endl;

		xmlNodeSetPtr productNodes = xmlXPathEvalExpression( (const xmlChar *) "./product", xpathCtx ) ->nodesetval;
		int numOfProducts = (productNodes) ? productNodes->nodeNr : 0;
		PROTECT(PRODUCTS = allocVector(STRSXP,numOfProducts));
		PROTECT(PSTOIC = allocVector(REALSXP,numOfProducts));

		for (int p = 0;p < numOfProducts;p++) {
			temp = get_attr(productNodes->nodeTab[p], "name");
			SET_STRING_ELT(PRODUCTS,p,mkChar(temp+4));
			REAL(PSTOIC)[p] = NA_REAL;
		}
		SET_VECTOR_ELT(REACTION,4,PRODUCTS); SET_STRING_ELT(REACTIONNAMES,4,mkChar("products"));
		SET_VECTOR_ELT(REACTION,5,PSTOIC);SET_STRING_ELT(REACTIONNAMES,5,mkChar("product.stoichiometry"));
		//cout << "Product level :|" << endl;

        SET_VECTOR_ELT(REACTION,6,R_NilValue);SET_STRING_ELT(REACTIONNAMES,6,mkChar("kinetics"));


        string geneXPath = ((string)"//entry[@type='gene' and @reaction='")+((string)name)+((string)"']");
        xmlNodeSetPtr geneNodes = xmlXPathEvalExpression(
        		(const xmlChar *) geneXPath.c_str(), xpathCtx ) ->nodesetval;
        int numOfModifiers = (geneNodes) ? geneNodes->nodeNr : 0;
		vector<string> genes;
		for (int m = 0;m < numOfModifiers;m++)
			genes = split( get_attr(geneNodes->nodeTab[m], "name"), ' ', genes);

		PROTECT(GENES = allocVector(STRSXP,genes.size()));
		PROTECT(KEGG_GENES = allocVector(STRSXP,genes.size()));
		PROTECT(NCBI_GENE = allocVector(STRSXP,genes.size()));
		for(size_t g=0; g<genes.size(); g++){
			SET_STRING_ELT(GENES,g, mkChar( genes[g].c_str() ));
			SET_STRING_ELT(KEGG_GENES,g, mkChar( genes[g].c_str() ));
			SET_STRING_ELT(NCBI_GENE,g, mkChar( genes[g].c_str()+4 ));
		}

		SET_VECTOR_ELT(REACTION,7,GENES); SET_STRING_ELT(REACTIONNAMES,7,mkChar("genes"));

		PROTECT(PATHWAY = allocVector(STRSXP,1));
		SET_STRING_ELT(PATHWAY,0,mkChar(pathwayTitle));
		SET_VECTOR_ELT(REACTION,8,PATHWAY); SET_STRING_ELT(REACTIONNAMES,8,mkChar("pathway"));

		// Set MIRIAM idnetifiers: kegg.pathway, kegg reaction, kegg.compound, kegg.genes, ncbi.gene
		PROTECT(KEGG_PATHWAY = allocVector(STRSXP,1));
		SET_STRING_ELT(KEGG_PATHWAY,0,mkChar(pathwayId));
		SET_VECTOR_ELT(REACTION,9,KEGG_PATHWAY); SET_STRING_ELT(REACTIONNAMES,9,mkChar("miriam.kegg.pathway"));

		std::vector<std::string> kegg_reaction = split(name, ' ');
		PROTECT(KEGG_REACTION = allocVector(STRSXP, kegg_reaction.size() ));
		for(size_t kr=0; kr<kegg_reaction.size(); kr++)
			SET_STRING_ELT(KEGG_REACTION,kr, mkChar(kegg_reaction[kr].c_str() +3));
		SET_VECTOR_ELT(REACTION,10,KEGG_REACTION); SET_STRING_ELT(REACTIONNAMES,10,mkChar("miriam.kegg.reaction"));

		SET_VECTOR_ELT(REACTION,11,KEGG_GENES); SET_STRING_ELT(REACTIONNAMES,11,mkChar("miriam.kegg.genes"));
		SET_VECTOR_ELT(REACTION,12,NCBI_GENE); SET_STRING_ELT(REACTIONNAMES,12,mkChar("miriam.ncbigene"));

		setAttrib(REACTION,R_NamesSymbol,REACTIONNAMES);
        SET_VECTOR_ELT(REACTIONLIST,i,REACTION);

        UNPROTECT(14);

    }

    setAttrib(REACTIONLIST,R_NamesSymbol,ID);
    UNPROTECT(2);

	/* Cleanup */
	xmlXPathFreeObject(nodes);
	xmlXPathFreeContext(xpathCtx);
	xmlFreeDoc(doc);

    //cout << "RacList returns :|" << endl;
    return(REACTIONLIST);
}

SEXP readkgml_sign(SEXP FILENAME, SEXP EXPAND_COMPLEXES, SEXP VERBOSE) {
	handle_segfault_KGML();

	bool expand_complexes = LOGICAL(EXPAND_COMPLEXES)[0];
	bool verbose = LOGICAL(VERBOSE)[0];

	vector<string> vertices;
	vector<int> edges;
	vector< vector<string> > attr;
	vector< vector<string> > pathway_attr;

	for(int i=0; i< LENGTH(FILENAME); i++){
		const char *filename = CHAR(STRING_ELT(FILENAME,i));
		readkgml_sign_int(filename, vertices, edges, attr, pathway_attr, expand_complexes, verbose);
	}

	SEXP VERTICES, V_NAMES,EDGES, E_ATTR;
	PROTECT(VERTICES = allocVector(VECSXP,vertices.size() ));
	PROTECT(V_NAMES = allocVector(STRSXP,vertices.size() ));
	PROTECT(EDGES = allocVector(INTSXP, edges.size() ));
	PROTECT(E_ATTR = allocVector(VECSXP, attr.size() ));

	/* Store vertex names and attributes */
	for(size_t vit = 0; vit < vertices.size(); vit++){
		SEXP V_ATTR, ATTR_NAMES;
		const char* v_name = vertices[vit].c_str();
		SET_STRING_ELT(V_NAMES, vit, mkChar( v_name ));

		// Compounds may be vertices in a signaling network via PCrel.
		if( strstr(v_name, "cpd") ){
			SEXP KEGG_CPD;
			PROTECT( V_ATTR = allocVector(VECSXP, 3) );
			PROTECT( ATTR_NAMES = allocVector(STRSXP, 3) );

			PROTECT(KEGG_CPD= allocVector(STRSXP,1));
			SET_STRING_ELT(KEGG_CPD, 0, mkChar( v_name +4 ));	SET_VECTOR_ELT(V_ATTR,0,KEGG_CPD);
			SET_STRING_ELT(ATTR_NAMES, 0, mkChar("miriam.kegg.compound"));

			UNPROTECT(1);//KEGG_CPD
		}else{
			SEXP KEGG_GENES, NCBI_GENE;
			PROTECT( V_ATTR = allocVector(VECSXP, 4) );
			PROTECT( ATTR_NAMES = allocVector(STRSXP, 4) );

			vector<string> genes = split( v_name, ' ');
			PROTECT(KEGG_GENES = allocVector(STRSXP,genes.size()));
			PROTECT(NCBI_GENE = allocVector(STRSXP,genes.size()));
			for(size_t g=0; g<genes.size(); g++){
				SET_STRING_ELT(KEGG_GENES,g, mkChar( genes[g].c_str() ));
				SET_STRING_ELT(NCBI_GENE,g, mkChar( genes[g].c_str() +4 ));
			}

			SET_VECTOR_ELT(V_ATTR,0,KEGG_GENES); SET_STRING_ELT(ATTR_NAMES,0,mkChar("miriam.kegg.genes"));
			SET_VECTOR_ELT(V_ATTR,1,NCBI_GENE); SET_STRING_ELT(ATTR_NAMES,1,mkChar("miriam.ncbigene"));
			UNPROTECT(2); //KEGG_GENES, NCBI_GENES
		}

		SEXP PATHWAY, KEGG_PATHWAY;
		PROTECT(PATHWAY = allocVector(STRSXP, pathway_attr[vit].size()/2));
		PROTECT(KEGG_PATHWAY = allocVector(STRSXP, pathway_attr[vit].size()/2));

		for(size_t p_at = 0; p_at<pathway_attr[vit].size()/2; p_at++){
			SET_STRING_ELT(KEGG_PATHWAY, p_at, mkChar( pathway_attr[vit][p_at*2].c_str() ));
			SET_STRING_ELT(PATHWAY, p_at, mkChar( pathway_attr[vit][p_at*2+1].c_str() ));
		}

		int pathway_pos = strstr(v_name, "cpd") ? 1 : 2;
		SET_VECTOR_ELT(V_ATTR, pathway_pos, KEGG_PATHWAY);
		SET_STRING_ELT(ATTR_NAMES, pathway_pos, mkChar("miriam.kegg.pathway"));
		SET_VECTOR_ELT(V_ATTR, pathway_pos+1, PATHWAY);
		SET_STRING_ELT(ATTR_NAMES, pathway_pos+1, mkChar("pathway"));
		UNPROTECT(2); //PATHWAY, KEGG_PATHWAY

		setAttrib(V_ATTR,R_NamesSymbol,ATTR_NAMES);
		SET_VECTOR_ELT(VERTICES,vit,V_ATTR);

		UNPROTECT(2); //V_ATTR, ATTR_NAMES
	}

	setAttrib(VERTICES,R_NamesSymbol,V_NAMES);

	/* Store edges and their attributes */
	for(size_t eit=0; eit<edges.size(); eit++)
		INTEGER(EDGES)[eit] = edges[eit]+1;

	for(size_t at_it=0; at_it<attr.size(); at_it++){
		SEXP EIT_ATTR, ATTR_NAMES;
		PROTECT( EIT_ATTR = allocVector(STRSXP, attr[at_it].size()) );
		PROTECT( ATTR_NAMES = allocVector(STRSXP, attr[at_it].size()) );

		for(size_t i=0; i<attr[at_it].size(); i++){
			const char* attr_i = attr[at_it][i].c_str();
			if(strstr(attr_i, "cpd")){
				SET_STRING_ELT(EIT_ATTR, i, mkChar( attr_i +4 )); SET_STRING_ELT(ATTR_NAMES, i, mkChar("miriam.kegg.compound"));
			}else{
				SET_STRING_ELT(EIT_ATTR, i, mkChar( attr_i )); SET_STRING_ELT(ATTR_NAMES, i, mkChar("type"));
			}
		}
		setAttrib(EIT_ATTR,R_NamesSymbol,ATTR_NAMES);
		SET_VECTOR_ELT(E_ATTR,at_it,EIT_ATTR);
		UNPROTECT(2);
	}

	SEXP RESULT;
	PROTECT( RESULT = allocVector(VECSXP, 3) );
	SET_VECTOR_ELT(RESULT,0,VERTICES);
	SET_VECTOR_ELT(RESULT,1,EDGES);
	SET_VECTOR_ELT(RESULT,2,E_ATTR);

	UNPROTECT(5);

	return(RESULT);
}

void readkgml_sign_int(const char* filename,
								vector<string> &vertices,
								vector<int> &edges,
								vector< vector<string> > &attr,
								vector< vector<string> > &pathway_attr,
								bool expand_complexes, bool verbose)
{
	xmlDocPtr doc;
	xmlXPathContextPtr xpathCtx;
	xmlXPathObjectPtr nodes;

	if(verbose)	Rprintf("Processing KGML file: %s",filename);

	/* Load XML document */
	doc = xmlParseFile(filename);
	if (doc == NULL) {
		Rf_warningcall(mkChar(filename), "Unable to parse file.");
		if(verbose)	Rprintf(": Error.\n");
		return;
	}

	//Check if the xml file has a KEGG DTD System.
	/* Check it is a kegg pathway file */
	if(doc->intSubset == NULL ||
	   strcmp( (char *) (doc->intSubset->name), "pathway") != 0 )
	   //strncmp( (char *) (doc->intSubset->SystemID), "http://www.kegg.jp/kegg/", 24) !=0)
	{
		Rf_warningcall(mkChar(filename), "File is not KEGG pathway file.");
		xmlFreeDoc(doc);
		if(verbose)	Rprintf(": Error.\n");
		return;
	}

	/* Get pathway information :*/
	xmlNodePtr pathway =  xmlDocGetRootElement(doc);
	if(!pathway || strcmp( (char *) (pathway->name), "pathway") != 0){
		Rf_warningcall(mkChar(filename), "No pathways in file.");
		xmlXPathFreeContext(xpathCtx);
		xmlFreeDoc(doc);
		if(verbose)	Rprintf(": Error.\n");
		return;
	}

	vector<string> pathway_info;
	const char* pathwayId = get_attr(pathway, "name");
	if(!pathwayId){
		Rf_warningcall(mkChar(filename), "Pathway ID not found in file. Using file name instead.");
		pathwayId = filename;
	}else{
		pathwayId +=5; //Remove "path:" leading characters//
	}

	const char* pathwayTitle = get_attr(pathway, "title");
	if(!pathwayTitle){
		Rf_warningcall(mkChar(pathwayId), "Pathway title not found in file.");
		pathwayTitle = "";
	}
	if(verbose)	Rprintf(" \"%s\"",pathwayTitle);

	/* Create xpath evaluation context */
	xpathCtx = xmlXPathNewContext(doc);
	if(xpathCtx == NULL) {
		Rf_warningcall(mkChar(filename), "Unable to create new XPath context.");
		xmlFreeDoc(doc);
		if(verbose)	Rprintf(": Error.\n");
		return;
	}

	/* Evaluate xpath expression */
	nodes = xmlXPathEvalExpression((xmlChar *) "//relation", xpathCtx);
	if(nodes == NULL || nodes->nodesetval == NULL || nodes->nodesetval->nodeNr == 0) {
		Rf_warningcall(mkChar(pathwayId), "Pathway contains no Protein-protein relationships.");
		xmlXPathFreeContext(xpathCtx);
		xmlFreeDoc(doc);
		if(verbose)	Rprintf(": Error.\n");
		return;
	}

	/* Parse XML Reactions*/
	xmlNodePtr curRelation;
    int size = (nodes->nodesetval) ? nodes->nodesetval->nodeNr : 0;

    if(verbose)	Rprintf(": %d gene relations found.\n",size);

    /* Looping over "relations" */
    for(int i = 0; i < size; ++i) {
		curRelation = nodes->nodesetval->nodeTab[i];
		char* type = get_attr(curRelation, "type");
		if(!type || strcmp(type, "maplink") == 0 )
			continue;

		// Get gene names for entry1 and entry2
		vector<string> p1,p2; //Holder objects for all gene names in this "relation"

		char* entry1 = get_attr(curRelation, "entry1");
		char* p1_name = entry1 ? attr_by_id(entry1, "name", xpathCtx) : NULL;
		if(p1_name && strcmp(p1_name, "undefined") == 0 ){
			p1_name = get_group_components(entry1, xpathCtx);
		}
		if(!p1_name) continue;

		char* entry2 = get_attr(curRelation, "entry2");
		char* p2_name = entry2 ? attr_by_id(entry2, "name", xpathCtx) : NULL;
		if(p2_name && strcmp(p2_name, "undefined") == 0 ){
			p2_name = get_group_components(entry2, xpathCtx);
		}
		if(!p2_name) continue;


		/* If complexes are expanded, each gene is a separate vertex.
		 * Otherwise, all genes participating in the relation are kept in a single vertex.
		 */
		if(expand_complexes){
			p1 = split(p1_name, ' ');
			p2 = split(p2_name, ' ');
		}else{
			p1.push_back(p1_name);
			p2.push_back(p2_name);
		}

		// Check if p1, p2 are already in our stack.
		vector<size_t> p1_pos, p2_pos;
		for(size_t j = 0; j < p1.size(); j++){
			p1_pos.push_back( elem_pos( vertices, p1[j] ) );

			if(p1_pos[j] == vertices.size())
				vertices.push_back(p1[j]);
		}

		for(size_t k = 0; k < p2.size(); k++){
			p2_pos.push_back( elem_pos( vertices, p2[k] ) );

			if(p2_pos[k] == vertices.size())
				vertices.push_back(p2[k]);
		}

		/* Setting pathway attributes for all added vertices */
		//making sure pathway and vertices vectors are of the same size//
		for(size_t p_attr = pathway_attr.size(); p_attr < vertices.size(); p_attr++)
			{ pathway_attr.push_back(vector<string>()); }

		//Adding this pathway as attribute, if it's not already added//
		string pid = pathwayId;
		for(size_t j=0; j<p1_pos.size(); j++){
			if( !elem_in_vector(pathway_attr[ p1_pos[j] ], pid ) ){
				pathway_attr[ p1_pos[j] ].push_back(pathwayId);
				pathway_attr[ p1_pos[j] ].push_back(pathwayTitle);
			}
		}

		for(size_t k=0; k<p2_pos.size(); k++){
			if( !elem_in_vector(pathway_attr[ p2_pos[k] ], pid ) ){
				pathway_attr[ p2_pos[k] ].push_back(pathwayId);
				pathway_attr[ p2_pos[k] ].push_back(pathwayTitle);
			}
		}

		/* Edges to connect the added vertices, and their attributes */
		// Relation parsing depennds on its type //
		if( strcmp(type, "PPrel") == 0 || strcmp(type, "GErel") == 0 || strcmp(type, "PCrel") == 0 ){
			// Add all combinatiosn from p1-> p2 as edges.
			for(size_t j=0; j<p1_pos.size(); j++){
				for(size_t k=0; k<p2_pos.size(); k++){
					edges.push_back(p1_pos[j]);	edges.push_back(p2_pos[k]);
				}
			}

			vector<string> e_attr;

			xpathCtx->node = curRelation;
			xmlNodeSetPtr subtype = xmlXPathEvalExpression( (const xmlChar *) "./subtype", xpathCtx ) ->nodesetval;
			int numOfattr = (subtype) ? subtype->nodeNr : 0;

			for (int a = 0;a < numOfattr;a++){
				xmlNodePtr sub_node = subtype->nodeTab[a];
				char* subtype_name = get_attr(sub_node, "name");

				if(!subtype_name)
					continue;

				if(strcmp(subtype_name, "compound") == 0){
					char* cpd_name = attr_by_id(get_attr(sub_node, "value"), "name", xpathCtx);
					e_attr.push_back(cpd_name);
				}else{
					e_attr.push_back(subtype_name);
				}
			}

			// Add the same attribute for all added edges.
			for(size_t l=0; l<p1_pos.size()*p2_pos.size(); l++)
				attr.push_back(e_attr);
		}
		else if(strcmp(type, "ECrel") == 0 ){
			/* ECrel indicated participation in 2 succesive reactions
			 * For ECrel, KGML deson't respect the direction of the relation
			 * Here, I will try to find whether it's entry1->entry2, or the reverse.
			 * Below, p1 particpates in r1, and p2 in r2, and the shared compound is cpd.
			 */
			char* cpd_id = get_attr(curRelation->children->next, "value");
			char* cpd = attr_by_id(cpd_id, "name", xpathCtx); //Cpd name
			if(!cpd) continue;
			char* r1_name = attr_by_id(entry1, "reaction",xpathCtx);
			xmlNodePtr r1node = r1_name ? node_by_attr_val("name", r1_name, "reaction",xpathCtx) : NULL;
			if(!r1node)continue;

			bool r1_rev = strcmp(get_attr(r1node, "type"), "reversible") == 0;
			bool r1_cpd = false; //If R1->Cpd (compound is a product of R1).

			if(!r1_rev){
				xpathCtx->node = r1node;
				string childXPath = ((string)"./*[@name='")+((string)cpd)+((string)"']");
				xmlNodeSetPtr children = xmlXPathEvalExpression(
							(const xmlChar *) childXPath.c_str(), xpathCtx ) ->nodesetval;
				char* role = children && children->nodeNr >0 ? (char*)children->nodeTab[0]-> name : NULL;
				if(!role) continue;

				if(strcmp( role , "product") == 0)
						r1_cpd = true;
				else{ r1_cpd = false;}
			}// !r1_rev

			char* r2_name = attr_by_id(entry2, "reaction",xpathCtx);
			xmlNodePtr r2node = r2_name ? node_by_attr_val("name", r2_name, "reaction",xpathCtx) : NULL;
			if(!r2node)continue;

			bool r2_rev = strcmp(get_attr(r2node, "type"), "reversible") == 0;
			bool r2_cpd = false;
			if(!r2_rev){
				xpathCtx->node = r2node;
				string childXPath = ((string)"./*[@name='")+((string)cpd)+((string)"']");
				xmlNodeSetPtr children = xmlXPathEvalExpression(
							(const xmlChar *) childXPath.c_str(), xpathCtx ) ->nodesetval;

				char* role = children && children->nodeNr >0 ? (char*)children->nodeTab[0]-> name : NULL;
				if(!role) continue;

				if(strcmp( role , "product") == 0)
						r2_cpd = true;
				else{ r2_cpd = false;}
			}// !r2_rev

			/* the order of r1, r2 is:
			 * r1 -> r2 if cpd is a product of r1 and substrate of 2
			 * meaning r1_cpd=true, r2_cpd=false.
			 * The opposite is also true.
			 * The value of r1_cpd doesn't matter if r1 reversible.
			 */
			if((r1_rev || r1_cpd) && (r2_rev || !r2_cpd)){
				for(size_t j=0; j<p1_pos.size(); j++){
					for(size_t k=0; k<p2_pos.size(); k++){
						edges.push_back(p1_pos[j]);	edges.push_back(p2_pos[k]);
						vector<string> e_attr;	e_attr.push_back(cpd);
						attr.push_back(e_attr);
					}
				}
			}// R1->R2

			if((r1_rev || !r1_cpd) && (r2_rev || r2_cpd)){
				for(size_t j=0; j<p1_pos.size(); j++){
					for(size_t k=0; k<p2_pos.size(); k++){
						edges.push_back(p2_pos[k]);	edges.push_back(p1_pos[j]);
						vector<string> e_attr;	e_attr.push_back(cpd);
						attr.push_back(e_attr);
					}
				}
			}// R2->R1 (order is reversed)
		}// End ECrel
    }// End for(relations)
}//kgml_sig_int

xmlNodePtr node_by_id(char* id, const char* type, xmlXPathContextPtr &xpathCtx){
	string XPath = ((string)"//")+ ((string)type)+ ((string)"[@id='")+ ((string)id)+ ((string)"']");
	xmlNodeSetPtr node_set = xmlXPathEvalExpression(
				(const xmlChar *) XPath.c_str(), xpathCtx ) ->nodesetval;
	return( node_set && node_set->nodeNr >0 ? node_set->nodeTab[0] : NULL);
}

xmlNodePtr node_by_attr_val(const char* attr, char* val, const char* type, xmlXPathContextPtr &xpathCtx){
	string XPath = ((string)"//")+ ((string)type)+ ((string)"[@")+
				((string)attr)+ ((string)"='")+ ((string)val)+ ((string)"']");
	xmlNodeSetPtr node_set = xmlXPathEvalExpression(
				(const xmlChar *) XPath.c_str(), xpathCtx ) ->nodesetval;
	return( node_set && node_set->nodeNr >0 ? node_set->nodeTab[0] : NULL);
}
char* get_attr(xmlNodePtr node, const char* attr_name){
	return( (char *) xmlGetProp(node,(const xmlChar *)attr_name) );
}
char* attr_by_id(char* id, const char* attr_name,xmlXPathContextPtr &xpathCtx){
	xmlNodePtr node = node_by_id(id, "entry", xpathCtx);
	return( node? get_attr(node, attr_name) : NULL );
}

char* get_group_components(char* id, xmlXPathContextPtr &xpathCtx){
	xpathCtx->node = node_by_id(id, "entry", xpathCtx);
	xmlNodeSetPtr comp = xmlXPathEvalExpression( (const xmlChar *) "./component", xpathCtx ) ->nodesetval;
	int numOfcomp = (comp) ? comp->nodeNr : 0;

	string group = "";
	for (int a = 0;a < numOfcomp;a++){
		xmlNodePtr sub_node = comp->nodeTab[a];
		char* comp_id = get_attr(sub_node, "id");
		char* comp_name = comp_id ? attr_by_id(comp_id, "name", xpathCtx) : NULL;

		if(comp_name)
			group = group + ((string)" ") + ((string) comp_name);

	}
	return(strdup(group.c_str()));
}

template <class T>
size_t elem_pos(vector<T> v, T &e){
	return( find(v.begin(), v.end(), e) - v.begin() );
}
template <class T>
bool elem_in_vector(vector<T> v, T &e){
	return( elem_pos(v,e) < v.size() );
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}
#endif
