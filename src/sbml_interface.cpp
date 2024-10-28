#ifdef HAVE_SBML
#include "sbml_interface.h"

template <class T>
void free_vec(vector<T> &v){
	vector<T>().swap(v);
}

SEXP readsbmlfile(SEXP FILENAME, SEXP ATTR_TERMS, SEXP VERBOSE) {
	handle_segfault_SBML();

	SEXP SPECIESFRAME, REACTIONLIST, OUT,NAMES;
	const char *filename = CHAR(STRING_ELT(FILENAME,0));

	vector<string> attr_terms;
	for(int i=0; i<LENGTH(ATTR_TERMS); i++)
	  attr_terms.push_back( CHAR(STRING_ELT(ATTR_TERMS,i)) );

	bool verbose = LOGICAL(VERBOSE)[0];

	SBMLDocument* document = readSBML(filename);
	unsigned int errors = document->getNumErrors();

	if(verbose){
		Rprintf("Processing SBML file: %s",filename);
		Rprintf( ", SBML level %d",document->getLevel());
		Rprintf( " version %d",document->getVersion());
	}

	Model *model = document->getModel();
	if(!model){
		if(verbose)	Rprintf(": Error.\n");
		Rf_warningcall(mkChar(filename), "No model in file");
		return(NEW_LIST(2));
	}

	if(errors>0){
		ostringstream message;
		for(unsigned e=0; e<errors; e++){
			const SBMLError *err = document->getError(e);
			message<<"line "<< err->getLine() <<": "<< err->getShortMessage() << "\n";
			if(err->getErrorId() == NotSchemaConformant ){
				if(verbose)	Rprintf(": Error.\n");
				Rf_warningcall(mkChar(filename), "%s",  message.str().c_str());
				return(R_NilValue);
			}
		}//loop over errors.
		Rf_warningcall(mkChar(filename), "%s", message.str().c_str());
	}




	vector<string> species;
	PROTECT( REACTIONLIST = getReactionList(model, attr_terms, species, verbose) );
	PROTECT( SPECIESFRAME = getSpeciesFrame(model, species, attr_terms) );


	//cout << "C++ returns :|" << endl;
	PROTECT( OUT = NEW_LIST(2) );
	PROTECT( NAMES = NEW_STRING(2) );
	SET_VECTOR_ELT(OUT,0,REACTIONLIST); SET_STRING_ELT(NAMES,0,mkChar("reactions"));
	SET_VECTOR_ELT(OUT,1,SPECIESFRAME); SET_STRING_ELT(NAMES,1,mkChar("species"));

	setAttrib(OUT,R_NamesSymbol,NAMES);

	UNPROTECT(4);

	return(OUT);
}

SEXP getReactionList(Model *model, const vector<string> &attr_terms, vector<string> &species, bool verbose) {
    ListOfReactions *reactions = model->getListOfReactions();
    ListOfSpecies *speciesList = model->getListOfSpecies();
    ListOfCompartments *compList = model->getListOfCompartments();
    
    if(verbose)	Rprintf(": %d reactions found.\n", reactions->size());
    
    SEXP REACTIONLIST,ID;
    PROTECT(REACTIONLIST = allocVector(VECSXP, reactions->size()));
    PROTECT(ID = NEW_STRING(reactions->size()));
    
    SEXP PATHWAY;
    PROTECT(PATHWAY = NEW_STRING(1));
    SET_STRING_ELT(PATHWAY,0, mkChar(model->getName().c_str()) );

    for (unsigned i = 0;i < reactions->size();i++) {
    	Reaction *ri = reactions->get(i);

    	// Attributes will include all associated with the reaction node and its modifiers
    	vector< vector<string> > attr;
    	vector<string> attr_names;
    	get_MIRIAM(ri->getAnnotation(), attr_terms, attr, attr_names);

    	// Compartment is inhireted from modifiers.
    	vector< vector<string> > comp_attr;
		vector<string> comp_attr_names, compartment, comp_name;
		vector<string> comp_attr_terms = attr_terms;
		comp_attr_terms.push_back("go");


		SET_STRING_ELT(ID,i,mkChar(ri->getId().c_str()));

        SEXP REACTION,REACTIONNAMES;
        SEXP NAME, REVERSIBLE, REACTANTS, RSTOIC, PRODUCTS, PSTOIC, GENES, KINETICS,KNAMES, COMPARTMENT, COMP_NAME;

        PROTECT( NAME = NEW_STRING(1) );
        SET_STRING_ELT(NAME,0, mkChar(ri->getName().c_str()) );
        
        PROTECT( REVERSIBLE = NEW_LOGICAL(1) );
        LOGICAL(REVERSIBLE)[0] = ri->getReversible();

        int numOfReactants = ri->getNumReactants();
        PROTECT( REACTANTS = NEW_STRING(numOfReactants) );
        PROTECT( RSTOIC = NEW_NUMERIC(numOfReactants) );
        for (int r = 0;r < numOfReactants;r++) {
        	const string sp = ri->getReactant(r)->getSpecies();
        	add_elem( species, sp );

            SET_STRING_ELT(REACTANTS,r,mkChar(sp.c_str()));
            REAL(RSTOIC)[r] = ri->getReactant(r)->getStoichiometry();
        }
        //cout << "Reactant level :|" << endl;
        
        int numOfProducts = ri->getNumProducts();
        PROTECT( PRODUCTS = NEW_STRING(numOfProducts) );
        PROTECT( PSTOIC = NEW_NUMERIC(numOfProducts) );
        for (int p = 0;p < numOfProducts;p++) {
        	const string sp = ri->getProduct(p)->getSpecies();
			add_elem( species, sp );

        	SET_STRING_ELT(PRODUCTS,p,mkChar(sp.c_str()));
            REAL(PSTOIC)[p] = ri->getProduct(p)->getStoichiometry();
        } 
        //cout << "Product level :|" << endl;

        //Kinetic law
        KineticLaw *kinetics = ri->getKineticLaw();
        int knum = kinetics ? kinetics->getNumParameters() : 0;
		//cout << "Kinetic for level :|" << knum << endl;
		PROTECT( KINETICS = NEW_LIST(knum) );
        PROTECT( KNAMES = NEW_STRING(knum) );

        for (int k = 0;k < knum;k++) {
           SEXP value;
           PROTECT( value = NEW_NUMERIC(1) );
           REAL(value)[0] = kinetics->getParameter(k)->getValue();
           SET_STRING_ELT(KNAMES,k,mkChar(kinetics->getParameter(k)->getId().c_str()));
           SET_VECTOR_ELT(KINETICS,k,value);
           UNPROTECT(1);
        }
        setAttrib(KINETICS,R_NamesSymbol,KNAMES);
        //cout << "kinetic level :|" << endl;

        int numOfModifiers = ri->getNumModifiers();
		PROTECT( GENES = NEW_STRING(numOfModifiers) );
		for (int m = 0;m < numOfModifiers;m++) {
			Species *sp = speciesList->get( ri->getModifier(m)->getSpecies() );
			get_MIRIAM(sp->getAnnotation(), attr_terms, attr, attr_names);

			//Compartment info and attributes.
			Compartment *comp = compList->get( sp->getCompartment() );
			get_MIRIAM(comp->getAnnotation(), comp_attr_terms, comp_attr, comp_attr_names);
			add_elem(compartment, comp->getId());	add_elem(comp_name, comp->getName());

			SET_STRING_ELT(GENES,m,mkChar(sp->getName().c_str()));
		}

		PROTECT(COMPARTMENT = NEW_STRING(compartment.size()));
		PROTECT(COMP_NAME = NEW_STRING(compartment.size()));
		for(size_t a=0; a< compartment.size(); a++){
			SET_STRING_ELT(COMPARTMENT,a, mkChar(compartment[a].c_str() ));
			SET_STRING_ELT(COMP_NAME,a, mkChar(comp_name[a].c_str() ));
		}

        PROTECT( REACTION = NEW_LIST(attr.size() + comp_attr.size() + 11) );
        PROTECT( REACTIONNAMES = NEW_STRING(attr.size() + comp_attr.size() + 11) );

        SET_VECTOR_ELT(REACTION,0,NAME); SET_STRING_ELT(REACTIONNAMES,0,mkChar("name"));
        SET_VECTOR_ELT(REACTION,1,REVERSIBLE); SET_STRING_ELT(REACTIONNAMES,1,mkChar("reversible"));
        SET_VECTOR_ELT(REACTION,2,REACTANTS); SET_STRING_ELT(REACTIONNAMES,2,mkChar("reactants"));
        SET_VECTOR_ELT(REACTION,3,RSTOIC);SET_STRING_ELT(REACTIONNAMES,3,mkChar("reactant.stoichiometry"));
        SET_VECTOR_ELT(REACTION,4,PRODUCTS); SET_STRING_ELT(REACTIONNAMES,4,mkChar("products"));
        SET_VECTOR_ELT(REACTION,5,PSTOIC);SET_STRING_ELT(REACTIONNAMES,5,mkChar("product.stoichiometry"));
        SET_VECTOR_ELT(REACTION,6,KINETICS);SET_STRING_ELT(REACTIONNAMES,6,mkChar("kinetics"));
		SET_VECTOR_ELT(REACTION,7,GENES); SET_STRING_ELT(REACTIONNAMES,7,mkChar("genes"));
		SET_VECTOR_ELT(REACTION,8,COMPARTMENT); SET_STRING_ELT(REACTIONNAMES,8,mkChar("compartment"));
		SET_VECTOR_ELT(REACTION,9,COMP_NAME); SET_STRING_ELT(REACTIONNAMES,9,mkChar("compartment.name"));
		SET_VECTOR_ELT(REACTION,10,PATHWAY); SET_STRING_ELT(REACTIONNAMES,10,mkChar("pathway"));

		// Put MIRIAM attributes in the list.
		int idx = 11;
		for(size_t a=0; a< attr.size(); a++){
			SEXP ATTR;
			PROTECT( ATTR = NEW_STRING(attr[a].size()) );
			for(size_t a_sub=0; a_sub< attr[a].size(); a_sub++)
				SET_STRING_ELT(ATTR, a_sub, mkChar( attr[a][a_sub].c_str() ));

			SET_VECTOR_ELT(REACTION, idx, ATTR);
			string at_name = ( ((string)"miriam.") + attr_names[a]);
			SET_STRING_ELT(REACTIONNAMES,idx,mkChar(at_name.c_str()));
			UNPROTECT(1);
			idx++;
		}
		//Add compartment attributes
		for(size_t a=0; a< comp_attr.size(); a++){
			SEXP ATTR;
			PROTECT( ATTR = NEW_STRING(comp_attr[a].size()) );
			for(size_t a_sub=0; a_sub< comp_attr[a].size(); a_sub++)
				SET_STRING_ELT(ATTR, a_sub, mkChar( comp_attr[a][a_sub].c_str() ));

			SET_VECTOR_ELT(REACTION, idx, ATTR);
			string at_name = ( ((string)"compartment.miriam.") + comp_attr_names[a]);
			SET_STRING_ELT(REACTIONNAMES,idx,mkChar( at_name.c_str()));
			UNPROTECT(1);
			idx++;
		}
		free_vec(attr); free_vec(attr_names);
		free_vec(comp_attr); free_vec(comp_attr_names); free_vec(comp_attr_terms);

        setAttrib(REACTION,R_NamesSymbol,REACTIONNAMES);
        SET_VECTOR_ELT(REACTIONLIST,i,REACTION);
        UNPROTECT(13);
    }
    
    setAttrib(REACTIONLIST,R_NamesSymbol,ID);
    UNPROTECT(3);
    //cout << "RacList returns :|" << endl;
    return(REACTIONLIST);
}

//
// Can be made more general if more info is present
//
SEXP getSpeciesFrame(Model *model, vector<string> species, const vector<string> &attr_terms) {
	SEXP SPECIESFRAME,ID;
	PROTECT( SPECIESFRAME = NEW_LIST(species.size()) );
	PROTECT( ID = NEW_STRING(species.size()) );

	for (size_t i = 0;i < species.size(); i++) {
		SET_STRING_ELT(ID,i, mkChar(species[i].c_str()));
		SET_VECTOR_ELT(SPECIESFRAME, i, get_species_info(model, species[i], attr_terms));
	}

  setAttrib(SPECIESFRAME,R_NamesSymbol,ID);

  UNPROTECT(2);
  return(SPECIESFRAME);
}

SEXP get_species_info(Model *model, const string species, const vector<string> &attr_terms){
	ListOfSpecies *speciesList = model->getListOfSpecies();
	ListOfCompartments *compList = model->getListOfCompartments();

	SEXP SP, SPNAMES;
	SEXP NAME,COMPARTMENT, COMP_NAME, PATHWAY;
	PROTECT(PATHWAY = NEW_STRING(1));
	SET_STRING_ELT(PATHWAY,0, mkChar(model->getName().c_str()) );

	Species *sp = speciesList->get( species );

	//Species attributes
	vector< vector<string> > attr;
	vector<string> attr_names;
	get_MIRIAM(sp->getAnnotation(), attr_terms, attr, attr_names);

	//Compartment info and attributes.
	Compartment *comp = compList->get( sp->getCompartment() );
	vector< vector<string> > comp_attr;
	vector<string> comp_attr_names;
	vector<string> comp_attr_terms = attr_terms;
	comp_attr_terms.push_back("go");
	get_MIRIAM(comp->getAnnotation(), comp_attr_terms, comp_attr, comp_attr_names);


	PROTECT(NAME = NEW_STRING(1));
	PROTECT(COMPARTMENT = NEW_STRING(1));
	PROTECT(COMP_NAME = NEW_STRING(1));
	SET_STRING_ELT(NAME,0, mkChar(sp->getName().c_str() ));
	SET_STRING_ELT(COMPARTMENT,0, mkChar(comp->getId().c_str() ));
	SET_STRING_ELT(COMP_NAME,0, mkChar(comp->getName().c_str() ));


	PROTECT( SP = NEW_LIST(attr.size() + comp_attr.size() + 4) );
	PROTECT( SPNAMES = NEW_STRING(attr.size() + comp_attr.size() + 4) );

	SET_VECTOR_ELT(SP, 0, NAME); SET_STRING_ELT(SPNAMES, 0, mkChar("name"));
	SET_VECTOR_ELT(SP, 1, COMPARTMENT); SET_STRING_ELT(SPNAMES, 1, mkChar("compartment"));
	SET_VECTOR_ELT(SP, 2, COMP_NAME); SET_STRING_ELT(SPNAMES, 2, mkChar("compartment.name"));
	SET_VECTOR_ELT(SP, 3, PATHWAY); SET_STRING_ELT(SPNAMES, 3, mkChar("pathway"));

	// Put MIRIAM attributes in the list.
	int idx = 4;
	for(size_t a=0; a< attr.size(); a++){
		SEXP ATTR;
		PROTECT( ATTR = NEW_STRING(attr[a].size()) );
		for(size_t a_sub=0; a_sub< attr[a].size(); a_sub++)
			SET_STRING_ELT(ATTR, a_sub, mkChar( attr[a][a_sub].c_str() ));

		SET_VECTOR_ELT(SP, idx, ATTR);
		string at_name = ( ((string)"miriam.") + attr_names[a]);
		SET_STRING_ELT(SPNAMES,idx,mkChar(at_name.c_str()));
		UNPROTECT(1);
		idx++;
	}
	//Add compartment attributes
	for(size_t a=0; a< comp_attr.size(); a++){
		SEXP ATTR;
		PROTECT( ATTR = NEW_STRING(comp_attr[a].size()) );
		for(size_t a_sub=0; a_sub< comp_attr[a].size(); a_sub++)
			SET_STRING_ELT(ATTR, a_sub, mkChar( comp_attr[a][a_sub].c_str() ));

		SET_VECTOR_ELT(SP, idx, ATTR);
		string at_name = ( ((string)"compartment.miriam.") + comp_attr_names[a]);
		SET_STRING_ELT(SPNAMES,idx,mkChar( at_name.c_str()));
		UNPROTECT(1);
		idx++;
	}
	free_vec(attr); free_vec(attr_names);
	free_vec(comp_attr); free_vec(comp_attr_names); free_vec(comp_attr_terms);

	//cout<<"species";
	setAttrib(SP,R_NamesSymbol,SPNAMES);
	UNPROTECT(6);
	return(SP);
}

SEXP readsbml_sign(SEXP FILENAME, SEXP ATTR_TERMS, SEXP VERBOSE){
	handle_segfault_SBML();

	vector<string> attr_terms;
	for(int i=0; i<LENGTH(ATTR_TERMS); i++)
	  attr_terms.push_back( CHAR(STRING_ELT(ATTR_TERMS,i)) );

	bool verbose = LOGICAL(VERBOSE)[0];

	vector<string> species;
	vector<size_t> edges, non_gene;
	vector<SEXP> info;

	for(int i=0; i<LENGTH(FILENAME); i++){
		const char *filename = CHAR(STRING_ELT(FILENAME,i));
		SBMLDocument* document = readSBML(filename);

		if(verbose){
			Rprintf("Processing SBML file: %s",filename);
			Rprintf( ", SBML level %d",document->getLevel());
			Rprintf( " version %d",document->getVersion());
		}

		Model *model = document->getModel();
		if(!model){
			if(verbose)	Rprintf(": Error.\n");
			Rf_warningcall(mkChar(filename), "No model in file");
			continue;
		}

		bool fatal=false;
		if(document->getNumErrors()>0){
			ostringstream message;
			for(unsigned e=0; e<document->getNumErrors(); e++){
				const SBMLError *err = document->getError(e);
				message<<"line "<< err->getLine() <<": "<< err->getShortMessage() << "\n";
				if(err->getErrorId() == NotSchemaConformant ){
					if(verbose)	Rprintf(": Error.\n");
					fatal=true; break;
				}
			}//loop over errors.
			Rf_warningcall(mkChar(filename), "%s", message.str().c_str());
		}
		if(fatal) continue;

		readsbml_sign_int(model, species, non_gene, info, edges, attr_terms, verbose);
	}//loop over fileList

	SEXP VERTICES, EDGES, ATTR, NONG, OUT, NAMES;
	PROTECT( VERTICES = NEW_STRING(species.size()) );
	PROTECT( EDGES = NEW_INTEGER(edges.size()) );
	PROTECT( ATTR = NEW_LIST(info.size()) );
	PROTECT( NONG = NEW_INTEGER(non_gene.size()) );

	for(size_t i=0; i<species.size(); i++)
		SET_STRING_ELT(VERTICES, i, mkChar( species[i].c_str() ));

	for(size_t i=0; i<edges.size(); i++)
		INTEGER(EDGES)[i] = edges[i]+1;

	for(size_t i=0; i<info.size(); i++)
		SET_VECTOR_ELT(ATTR, i, info[i] );

	for(size_t i=0; i<non_gene.size(); i++)
		INTEGER(NONG)[i] = non_gene[i]+1;

	PROTECT( OUT = NEW_LIST(4));
	PROTECT( NAMES = NEW_STRING(4));
	SET_VECTOR_ELT(OUT, 0, VERTICES); SET_STRING_ELT(NAMES, 0, mkChar("vertices"));
	SET_VECTOR_ELT(OUT, 1, EDGES);	SET_STRING_ELT(NAMES, 1, mkChar("edges"));
	SET_VECTOR_ELT(OUT, 2, ATTR);	SET_STRING_ELT(NAMES, 2, mkChar("attr"));
	SET_VECTOR_ELT(OUT, 3, NONG);	SET_STRING_ELT(NAMES, 3, mkChar("non.gene"));

	setAttrib(OUT,R_NamesSymbol,NAMES);
	UNPROTECT(info.size());
	UNPROTECT(6);
	return(OUT);
}

void readsbml_sign_int(Model *model, vector<string> &species, vector<size_t> &non_gene,
						vector<SEXP> &info, vector<size_t> &edges,
						const vector<string> &attr_terms, bool verbose)
{
	ListOfReactions *reactions = model->getListOfReactions();

	if(verbose)	Rprintf(": %d reactions found.\n", reactions->size());

	for (unsigned i = 0;i < reactions->size();i++) {
		Reaction *ri = reactions->get(i);

		vector<size_t> reactants, products, modifiers;
		for (unsigned r = 0;r < ri->getNumReactants();r++){
			string sp = ri->getReactant(r)->getSpecies();

			size_t pos = add_elem(species, sp);
			reactants.push_back(pos);

			if(pos==info.size()){
				SEXP INFO;
				PROTECT(INFO = get_species_info(model, sp, attr_terms));
				info.push_back(INFO);
			}
		}

		for (unsigned p = 0;p < ri->getNumProducts();p++){
			string sp = ri->getProduct(p)->getSpecies();

			size_t pos = add_elem(species, sp);
			products.push_back(pos);

			if(pos==info.size()){
				SEXP INFO;
				PROTECT(INFO = get_species_info(model, sp, attr_terms));
				info.push_back(INFO);
			}
		}


		for (unsigned m = 0;m < ri->getNumModifiers();m++) {
			string sp = ri->getModifier(m)->getSpecies();

			size_t pos = add_elem(species, sp);
			modifiers.push_back(pos);

			if(pos==info.size()){
				SEXP INFO;
				PROTECT(INFO = get_species_info(model, sp, attr_terms));
				info.push_back(INFO);
			}
		}

		if(ri->getNumModifiers()==0){
			size_t pos = add_elem(species, ri->getId());
			modifiers.push_back(pos);
			add_elem(non_gene, pos);

			SEXP INFO, NAME, NAME_;
			PROTECT( INFO = NEW_LIST(1) ); PROTECT( NAME = NEW_STRING(1) ); PROTECT( NAME_ = NEW_STRING(1) );
			SET_STRING_ELT(NAME,0, mkChar(ri->getName().c_str()) );
			SET_STRING_ELT(NAME_,0, mkChar("name") );
			SET_VECTOR_ELT(INFO, 0, NAME); SET_NAMES(INFO, NAME_);
			info.push_back(INFO);
			UNPROTECT(2);
		}

		/* Add edges from r->m->p */
		for(size_t m=0; m<modifiers.size(); m++){
			for(size_t r=0; r<reactants.size(); r++){
				edges.push_back(reactants[r]); edges.push_back(modifiers[m]);
			}//connect r->m
			for(size_t p=0; p<products.size(); p++){
				edges.push_back(modifiers[m]); edges.push_back(products[p]);
			}//connect m->p
		}

	}//loop over reactions
}

const char* URL_decode(const char* URL){
	SEXP x = EVAL(lang2(install("URLdecode"), mkString(URL)));
	return(CHAR(STRING_ELT(x,0)));
}

void get_MIRIAM(XMLNode* rdf, const vector<string> &terms, vector< vector<string> > &values, vector<string> &names){
	if(terms[0]=="none")
		return;

	string URI;
	List cvl = List();
	RDFAnnotationParser::parseRDFAnnotation(rdf, (List *) &cvl);

	for(unsigned i=0; i<cvl.getSize();i++){
		CVTerm *cv = (CVTerm *) cvl.get(i);

		// Extract attributes for biological qualifiers of types bqb:is, bqb:hasPart only.
		if(cv->getQualifierType() == BIOLOGICAL_QUALIFIER
				&& cv->getBiologicalQualifierType() != BQB_IS && cv->getBiologicalQualifierType() != BQB_HAS_PART)
			continue;

		for(size_t n=0;n< cv -> getNumResources();n++){
			URI= cv ->getResourceURI(n);
			for(size_t t=0; t<terms.size(); t++){
				string term = terms[t];
				int pos;

				if(terms[t]=="all"){
					pos = URI.find("identifiers.org")+16;
					if(pos<16)pos = URI.find("miriam")+7;
					if(pos<7) break;

					int term_end = find_if(URI.begin()+pos, URI.end(), not_alnum) - URI.begin();
					term = URI.substr(pos, term_end-pos);
				}else{
					pos = URI.find(term);
				}

				if(pos > 0){
					size_t term_pos = elem_pos(names, term);
					if( term_pos==names.size() ){
						names.push_back( term );
						values.push_back(vector<string>());
					}
					values[ term_pos ].push_back( URL_decode( URI.substr(pos + term.length() +1).c_str() ));
					break;
				}// If term is found
			}// loop over terms
		}//loop over URIs
	}// loop over CVTerms
}

template <class T>
size_t elem_pos(vector<T> v, const T &e){
	return( find(v.begin(), v.end(), e) - v.begin() );
}

template <class T>
size_t add_elem(vector<T> &v,const T &e){
	size_t pos = elem_pos(v,e);
	if(pos == v.size())
		v.push_back(e);

	return pos;
}
bool not_alnum(char c){
	return(!( isalnum(c) || c=='.'));
}
#endif

