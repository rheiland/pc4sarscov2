
int check_custom_data()
{
    static int nV_external = microenvironment.find_density_index( "virion" ); 
    static int nA_external = microenvironment.find_density_index( "assembled virion" ); 
    static int nINF1 = microenvironment.find_density_index( "interferon 1" );

    static int nV_internal = pCell->custom_data.find_variable_index( "virion" ); 
    static int nA_internal = pCell->custom_data.find_variable_index( "assembled_virion" ); 

	double additional_death_rate = pCell->custom_data["max_infected_apoptosis_rate"] ; 

}
