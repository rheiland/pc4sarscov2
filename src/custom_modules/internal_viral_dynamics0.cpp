#include "./internal_viral_dynamics.h" 

using namespace PhysiCell; 

std::string internal_virus_replication_version = "0.3.0"; 

Submodel_Information internal_viral_dynamics_info; 

double viral_RNA_min= 1.e6;
double viral_RNA_max= -1.e6;
double viral_protein_min= 1.e6;
double viral_protein_max= -1.e6;
double dR_min = 1.e6;
double dR_max = -1.e6;
double dA_min = 1.e6;
double dA_max = -1.e6;

void internal_virus_model_setup( void )
{
		// set version
	internal_viral_dynamics_info.name = "internal viral replication dynamics"; 
	internal_viral_dynamics_info.version = internal_virus_replication_version; 
		// set functions 
	internal_viral_dynamics_info.main_function = NULL; 
	internal_viral_dynamics_info.phenotype_function = internal_virus_model; 
	internal_viral_dynamics_info.mechanics_function = NULL; 
		// what microenvironment variables do I need? 

		// what custom data do I need? 
	internal_viral_dynamics_info.microenvironment_variables.push_back( "assembled virion" ); 	

	internal_viral_dynamics_info.cell_variables.push_back( "virion" ); // adhered, in process of endocytosis 
	internal_viral_dynamics_info.cell_variables.push_back( "uncoated_virion" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "viral_RNA" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "viral_protein" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "assembled_virion" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "export_virion" ); 

	internal_viral_dynamics_info.cell_variables.push_back( "virion_uncoating_rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "uncoated_to_RNA_rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "protein_synthesis_rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "virion_assembly_rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "virion_export_rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "max_RNA_replication_rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "RNA_replication_half" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "basal_RNA_degradation_rate" ); 

	// submodel_registry.register_model( internal_viral_dynamics_info ); 
	internal_viral_dynamics_info.register_model();
	
	return; 
}

void internal_virus_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	// bookkeeping -- find microenvironment variables we need

	static int nV_external = microenvironment.find_density_index( "virion" ); 
	static int nA_external = microenvironment.find_density_index( "assembled virion" ); 
	
	static int nV_internal = pCell->custom_data.find_variable_index( "virion" ); 
	static int nA_internal = pCell->custom_data.find_variable_index( "assembled_virion" ); 

	static int nUV = pCell->custom_data.find_variable_index( "uncoated_virion" ); 
	static int nR  = pCell->custom_data.find_variable_index( "viral_RNA" ); 
	static int nP  = pCell->custom_data.find_variable_index( "viral_protein" ); 	
	static int eP  = pCell->custom_data.find_variable_index( "export_virion" ); 

/*	
	static bool done = false; 
	extern Cell* pInfected; 
	if( pCell == pInfected && 1 == 0 )
	{
		std::cout << std::endl << "viral dynamics : " << __LINE__ << " " 
			<< phenotype.molecular.internalized_total_substrates[ nV_external ] << " " 
			<< phenotype.molecular.internalized_total_substrates[ nA_external ] << " " 
			<< pCell->custom_data[nV_internal] << " " 
			<< pCell->custom_data[nUV] << " " 
			<< pCell->custom_data[nR] << " " 
			<< pCell->custom_data[nP] << " " 	
			<< pCell->custom_data[nA_internal] << " " 
			<< std::endl; 		
	}
*/	
	
	// copy virions from "internalized variables" to "custom variables"
/*	
	pCell->custom_data[nV_internal] = 
		phenotype.molecular.internalized_total_substrates[nV_external]; 
	// this transfer is now handled in receptor dynamics 
*/		

  // This resets the internal assembled virion count 
  // so we are commenting it out 
	// pCell->custom_data[nA_internal] = 
	// 	phenotype.molecular.internalized_total_substrates[nA_external]; 
		
	// actual model goes here 

	// uncoat endocytosed virus
	double dV = dt * pCell->custom_data["virion_uncoating_rate"] * pCell->custom_data[nV_internal] ;
	if( dV > pCell->custom_data[nV_internal] )
	{ dV = pCell->custom_data[nV_internal]; } 
	pCell->custom_data[nV_internal] -= dV; 
	pCell->custom_data[nUV] += dV; 

	// convert uncoated virus to usable mRNA 
	double dR = dt * pCell->custom_data["uncoated_to_RNA_rate"] * pCell->custom_data[nUV]; 

    // gotta remove this from uncoated virions now befoe we add the replication
	if( dR > pCell->custom_data[nUV] )
	{ dR = pCell->custom_data[nUV]; }
    if (dR < dR_min)
    {
        dR_min = dR;
    }
    if (dR > dR_max)
    {
        dR_max = dR;
    }

	pCell->custom_data[nUV] -= dR; 
    // RNA replication post uncoated to RNA calc
	dR += dt * pCell->custom_data["max_RNA_replication_rate"] * pCell->custom_data[nR] /
              (pCell->custom_data[nR] + pCell->custom_data["RNA_replication_half"]);
    // RNA degradation
    dR -= dt * pCell->custom_data["basal_RNA_degradation_rate"] * pCell->custom_data[nR];

    // if( dR < -1*pCell->custom_data[nR] )
	// { dR = -1*pCell->custom_data[nR]; }

    // if (dR < 0.0) 
    // {
    //     // std::cout << "------- prevent dR < 0.0\n";
    //     dR = 0.0;  //rwh
    // }

	pCell->custom_data[nR] += dR; 
	// if (pCell->custom_data[nR] < 0.0)
    // {
	//     pCell->custom_data[nR] = 0.0;
    // }
    if (pCell->custom_data[nR] < viral_RNA_min)
    {
        viral_RNA_min = pCell->custom_data[nR];
    }
    if (pCell->custom_data[nR] > viral_RNA_max)
    {
        viral_RNA_max = pCell->custom_data[nR];
    }
	
	// use mRNA to create viral protein 
	double dP = dt * pCell->custom_data["protein_synthesis_rate"] * pCell->custom_data[nR];
	pCell->custom_data[nP] += dP; 
    // rwh: protein_synthesis_rate = 0.01
    if (pCell->custom_data[nP] < viral_protein_min)
    {
        viral_protein_min = pCell->custom_data[nP];
    }
    if (pCell->custom_data[nP] > viral_protein_max)
    {
        viral_protein_max = pCell->custom_data[nP];
    }

	// degrade protein 
	
	// assemble virus:  (rwh: nP = viral_protein) 
	double dA = dt * pCell->custom_data["virion_assembly_rate"] * pCell->custom_data[nP]; 
	// if( dA > pCell->custom_data[nP] )
	// { dA = pCell->custom_data[nP]; } 
    if (dA < dA_min)
    {
        dA_min = dA;
    }
    if (dA > dA_max)
    {
        dA_max = dA;
    }
	pCell->custom_data[nP] -= dA; 
	pCell->custom_data[nA_internal] += dA; 

    // if (int(PhysiCell_globals.current_time) % 100 == 0)
    // if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
    // if (PhysiCell_globals.current_time > 1440)
    // {
    //     std::cout << "-----" << __FUNCTION__ << ":  v_a_rate= " << pCell->custom_data["virion_assembly_rate"] << ", dA= " << dA << std::endl;
    //     std::cout << "            " << pCell->custom_data[nA_internal] << std::endl;
    // }
	
	// set export rate 
/*	
	phenotype.secretion.net_export_rates[nA_external] = 
		pCell->custom_data["virion_export_rate" ] * pCell->custom_data[nA_internal]; 
 
	// copy data from custom variables to "internalized variables" 
		
	phenotype.molecular.internalized_total_substrates[nV_external] = 
		pCell->custom_data[nV_internal];		
	phenotype.molecular.internalized_total_substrates[nA_external] = 
		pCell->custom_data[nA_internal];	
*/			
	double deP = dt * pCell->custom_data["virion_export_rate" ] * pCell->custom_data[nA_internal]; 
	if( deP > pCell->custom_data[nA_internal] )
	{ deP = pCell->custom_data[nA_internal]; } 
	pCell->custom_data[eP] += deP; 
	pCell->custom_data[nA_internal] -= deP; 
	//test int export
	
	double alpha1 = floor(pCell->custom_data[eP]);
	#pragma omp critical
	{ pCell->nearest_density_vector()[nV_external] += alpha1 / microenvironment.mesh.dV; }
	pCell->custom_data[eP] -= alpha1; 
	
	return; 
}

void internal_virus_model_symbolic( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int nV_external = microenvironment.find_density_index( "virion" ); 
	static int nA_external = microenvironment.find_density_index( "assembled virion" ); 
	
	static int nV_internal = pCell->custom_data.find_variable_index( "virion" ); 
	static int nA_internal = pCell->custom_data.find_variable_index( "assembled_virion" ); 

	static int nUV = pCell->custom_data.find_variable_index( "uncoated_virion" ); 
	static int nR  = pCell->custom_data.find_variable_index( "viral_RNA" ); 
	static int nP  = pCell->custom_data.find_variable_index( "viral_protein" ); 	
	static int eP  = pCell->custom_data.find_variable_index( "export_virion" ); 

  // This resets the internal assembled virion count so we are commenting it out 
	// pCell->custom_data[nA_internal] = 
	// 	phenotype.molecular.internalized_total_substrates[nA_external]; 
		
/*
	// ---- actual model goes here 
	// uncoat endocytosed virus
	double dV = dt * pCell->custom_data["virion_uncoating_rate"] * pCell->custom_data[nV_internal] ;
	if( dV > pCell->custom_data[nV_internal] )
	{ dV = pCell->custom_data[nV_internal]; } 

	pCell->custom_data[nV_internal] -= dV; 
	pCell->custom_data[nUV] += dV; 

	// convert uncoated virus to usable mRNA 
	double dR = dt * cell["uncoated_to_RNA_rate"] * cell[uncoated_virion]; 

    // gotta remove this from uncoated virions now befoe we add the replication
	if( dR > pCell->custom_data[nUV] )
	{ dR = pCell->custom_data[nUV]; }
    if (dR < dR_min)
    {
        dR_min = dR;
    }
    if (dR > dR_max)
    {
        dR_max = dR;
    }

	pCell->custom_data[nUV] -= dR; 
    // RNA replication post uncoated to RNA calc
	dR += dt * pCell->custom_data["max_RNA_replication_rate"] * pCell->custom_data[nR] /
              (pCell->custom_data[nR] + pCell->custom_data["RNA_replication_half"]);
    // RNA degradation
    dR -= dt * pCell->custom_data["basal_RNA_degradation_rate"] * pCell->custom_data[nR];

    // if( dR < -1*pCell->custom_data[nR] )
	// { dR = -1*pCell->custom_data[nR]; }

    // if (dR < 0.0) 
    // {
    //     // std::cout << "------- prevent dR < 0.0\n";
    //     dR = 0.0;  //rwh
    // }

	pCell->custom_data[nR] += dR; 
	// if (pCell->custom_data[nR] < 0.0)
    // {
	//     pCell->custom_data[nR] = 0.0;
    // }
    if (pCell->custom_data[nR] < viral_RNA_min)
    {
        viral_RNA_min = pCell->custom_data[nR];
    }
    if (pCell->custom_data[nR] > viral_RNA_max)
    {
        viral_RNA_max = pCell->custom_data[nR];
    }
	
	// use mRNA to create viral protein 
	double dP = dt * pCell->custom_data["protein_synthesis_rate"] * pCell->custom_data[nR];
	pCell->custom_data[viral_protein] += dP; 
    // rwh: protein_synthesis_rate = 0.01

	// degrade protein 
	
	// assemble virus:  (rwh: nP = viral_protein) 
	double dA = dt * cell[virion_assembly_rate] * cell[viral_protein]; 
	double dA = dt * cell[virion_assembly_rate] * cell[viral_protein]; 
	cell[viral_protein] -= dA; 
	cell[nA_internal] += dA; 

	// set export rate 
	double deP = dt * cell["virion_export_rate" ] * cell[nA_internal]; 
	if( deP > pCell->custom_data[nA_internal] )
	{ deP = pCell->custom_data[nA_internal]; } 
	pCell->custom_data[eP] += deP; 
	pCell->custom_data[nA_internal] -= deP; 
	//test int export
	
	double alpha1 = floor(pCell->custom_data[eP]);
	#pragma omp critical
	{ pCell->nearest_density_vector()[nV_external] += alpha1 / microenvironment.mesh.dV; }
	pCell->custom_data[eP] -= alpha1; 
	
*/			
	return; 
}
