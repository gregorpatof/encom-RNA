#include "STeM.h"


int main(int argc, char *argv[]) {

	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
	int all_t; /*Nombre d'atomes dans pdb*/
	int atom_t; /*Nombre de carbone CA*/
 	int help_flag = 1;
 	char file_name[500];
 	char check_name[500];
 	char eigen_name[100] = "eigen.dat";
 	int verbose = 0;
	int super_node = 1;
	int i;
	int nbr_mode = 2;
	int mode = 7;
	int lig = 0;
	int nconn;
	float ligalign = 0; // Flag/valeur pour aligner seulement les résidus dans un cutoff du ligand, 0, one le fait pas... > 0... le cutoff
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-sp",argv[i]) == 0) {super_node = 3;}
 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;} 
 		if (strcmp("-t",argv[i]) == 0) {strcpy(check_name,argv[i+1]);help_flag = 0;}
 		if (strcmp("-m",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);mode = temp;}
 		if (strcmp("-nm",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);nbr_mode = temp;}
 		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
 		if (strcmp("-ligc",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);ligalign = temp;}
 	}
 	
 	if (help_flag == 1) {
 		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-v\tVerbose\n-w\tWeight Vector\n-t\tInitial Value of template (negative value for random)\n\tIf Load Template, multiply the template\n-lt\tLoad input template\n-sp\tSuper Node Mode (CA, N, C)\n-kt\tPoid de l'angle entre les nodes (1)\n-kr\tPoid de la distance entre les nodes (1)\n-f\tFile to fit\n****************************\n");
 		return(0); 
 	} 

 	//***************************************************
 	//*													*
 	//*Build a structure contaning information on the pdb
 	//*													*
 	//***************************************************
 	
 	// Première strucutre
 	
 	all = count_atom(file_name);
 	nconn = count_connect(file_name);
 	if (verbose == 1) {printf("First file:%s\n",file_name);}
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_h[i]=(int *)malloc(4*sizeof(int));}
    
    assign_connect(file_name,connect_h);

	// Assign tous les atoms
	
	struct pdb_atom strc_all[all];
	atom = build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	
	if (verbose == 1) {printf("	Atom:%d\n",all);}
	check_lig(strc_all,connect_h,nconn,all);
	
	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}

	// Free Connect
		
	//for(i=0;i<nconn;i++) {printf("I:%d\n",i);free(connect_h[i]);}
	//free(connect_h);
	
	//Construit la structure a comparer
 	nconn = 0;
 	all_t = count_atom(check_name);
 	nconn = count_connect(check_name);
 	if (verbose == 1) {printf("Sec file:%s\n",check_name);}
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_t=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_t[i]=(int *)malloc(4*sizeof(int));}
    
    assign_connect(check_name,connect_t);
	
	// Assign tous les atoms
	
	struct pdb_atom strc_all_t[all_t];
	atom_t = build_all_strc(check_name,strc_all_t); // Retourne le nombre de Node
	if (verbose == 1) {printf("	Atom:%d\n",all_t);}
	check_lig(strc_all_t,connect_t,nconn,all_t);
	
	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom_t);}
	struct pdb_atom strc_node_t[atom_t];

	atom_t = build_cord_CA(strc_all_t, strc_node_t,all_t,lig,connect_t,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom_t);}
 
 	int align[atom];
 	int score = node_align(strc_node,atom,strc_node_t,atom_t,align);
 	if (verbose == 1) {printf("Score: %d/%d\n",score,atom);}
	//if (atom_t != atom) {printf("Not the same number of Node... Terminating\n");return(0);}
	
	if (score/atom < 0.8) {
 		printf("Low Score... Will try an homemade alignement !!!\n");
 		score = node_align_low(strc_node,atom,strc_node_t,atom_t,align);
 		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
 	}
		
	if (ligalign > 0) {

		score = node_align_lig(strc_node,atom,strc_node_t,atom_t,align,strc_all,all,strc_all_t,all_t,ligalign);
		
		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
	}
	//***************************************************
 	//*													*
 	//* Load eigenvector et eigenvalue					*
 	//*													*
 	//***************************************************
 	
 	if (verbose == 1) {printf("Loading Eigenvector\n");}
 	
 	gsl_vector *eval = gsl_vector_alloc(3*atom+10);
	gsl_matrix *evec = gsl_matrix_alloc (3*atom+10,3*atom+10);
 	
	
	load_eigen(eval,evec,eigen_name,3*atom);
	
	
	// Fit Using Eigenvalue
	if (nbr_mode == -1) {nbr_mode = 3*atom-7;}
	//fit(strc_node,strc_node_t,atom,all,strc_all,strc_all_t,evec,align);
	// nrg_rmsd(strc_node,atom,evec, align, nbr_mode, mode,eval);
	fit_vince(strc_node,strc_node_t,atom,all,strc_all,strc_all_t,evec,align,nbr_mode,mode,eval);
	//fit_math(strc_node,strc_node_t,atom,all,strc_all,strc_all_t,evec,align);
	
	//write_strc("target.pdb",strc_all_t,all);
 	//write_strc("final.pdb",strc_all,all);
	
    
    free(connect_h);
	free(connect_t);
	return(0);
 }
