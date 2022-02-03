Be sure PhysiCell_settings.xml "default" cell_def includes all 4 secretion params for *each* substrate.
e.g.,
                                        <substrate name="pro-pyroptosis cytokine">
                                                <secretion_target units="dimensionless substrate concentration">1</secretion_target>
                                                <uptake_rate units="1/min">0.0</uptake_rate>
                                                <secretion_rate units="1/min">0</secretion_rate>
                        <net_export_rate units="total substrate/min">0</net_export_rate>
                                        </substrate>

Be sure to change output folder:
	<omp_num_threads>4</omp_num_threads>
and:
        <folder>.</folder>		

Be sure to edit output intervals as you want (equal for now).

Be sure to include:
        <virtual_wall_at_domain_edge>true</virtual_wall_at_domain_edge>		

and:
    <initial_conditions>
		<cell_positions type="csv" enabled="false">
			<folder>./data</folder>
			<filename>cells.csv</filename>
		</cell_positions>
	</initial_conditions>		


If in the future there are <Dirichlet_boundary_condition enabled="true"..., be sure to specify the separate boundaries! 


--------------
~/git/pc4sarscov2/data$ python flatten_cell_def_xml.py PhysiCell_settings.xml
--> flat.xml, rename to whatever is desired.

Finally:
- remove "default" cell_def 
- remove "parent" attrib
- renumber all IDs for cell_defs, starting with 0 (really necessary??)

-->
~/git/pc4sarscov2/data$ grep "<cell_def" PhysiCell_settings_sarscov2.xml
	<cell_definitions>
		<cell_definition name="lung epithelium" ID="0" >
		<cell_definition name="CD8 Tcell" ID="1" >
		<cell_definition name="macrophage" ID="2" >
		<cell_definition name="neutrophil" ID="3" >
		<cell_definition name="DC" ID="4" >
		<cell_definition name="CD4 Tcell" ID="5" >
		<cell_definition name="fibroblast" ID="6" >
		<cell_definition name="residual" ID="7" >
~/git/pc4sarscov2/data$ 



