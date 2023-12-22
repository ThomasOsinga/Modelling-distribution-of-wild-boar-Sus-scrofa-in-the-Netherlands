


model Batchtest


import "../models/Boar_sounder_model_v_26-9.gaml"



global {
	
	int end_cycle <- 500;
		reflex save_raster when: cycle = 200 {
		ask suitability_cell {
			self.grid_value <- self.boar_time;
		}

		string file_name <- "../results/test/boartimeraster" + "_seed" + seed + ".tif";
		save suitability_cell to: file_name format: geotiff crs: 32631;
		do pause;
	}
	bool stop_sim { return (end_cycle = 500); } 
}


experiment batch_abstract type:batch virtual:true until:(time > end_cycle) {
	
}


experiment 'Run_5_simulations' parent: batch_abstract type: batch repeat: 5 parallel: 5 keep_seed: false until: world.stop_sim() or (time > end_cycle){
	
	// the reflex will be activated at the end of each run; in this experiment a run consists of the execution of 5 simulations (repeat: 5)
	reflex end_of_runs
	{
		ask 5 among simulations
		{
			string file_path <- "../Results/boar_sound_" +"seed_" + simulation.seed + ".csv";
			save (nb_male_adults + nb_female_adults + nb_male_yearlings + nb_female_yearlings + nb_piglets) to: file_path format:"csv";
		}
	}
}
