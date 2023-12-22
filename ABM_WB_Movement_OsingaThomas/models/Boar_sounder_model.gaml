
/**
* Name: AgentBasedBoarDistributionModel
* Author: Thomas Osinga
* Tags: ABM, Boar, Distribution, SDM, ASF
*/
model AgentBasedBoarDistributionModel

global {
	float seed <- 0.0; // setting the seed to 0.0 to randomize it for each run
 reflex save_raster when: cycle = 3000 {
	// This code block is creating an output raster to create 																				
 ask suitability_cell { // the utilization density - I replace grid_value with boar_time 

			self.grid_value <- self.boar_time; // because the grid value it exports is always grid_value
 }

		string file_name <- "../results/boartimeraster" + "_seed" + seed + ".tif";
		save suitability_cell to: file_name format: geotiff crs: 32631;
		do pause;
	}

	/* the following code blocks counts the number of animals each cycle to show in the graph
 * There is a difference between the first 4, as they are microspecies and thus have to be 'collected' from the sounder information
 * The other species can be directly counted by looking at the length of agent list
 */
	int nb_female_adults;
	int nb_male_yearlings;
	int nb_female_yearlings;
	int nb_piglets;
	int nb_sounders -> {length(sounder)};
	int nb_dispersing_females -> {length(dispersing_females)};
	int nb_dispersing_males -> {length(dispersing_males)};
	int nb_male_adults -> {length(adult_male)};

	action update_counts {
		nb_female_adults <- sum(sounder collect each.local_female_adults);
		nb_male_yearlings <- sum(sounder collect each.local_male_yearlings);
		nb_female_yearlings <- sum(sounder collect each.local_female_yearlings);
		nb_piglets <- sum(sounder collect each.local_piglets);
	}

	reflex {
		do update_counts;
	}

	/* here starts the initialization of the model parameters, these are the global parameters setting the scene for the model 
 * We start with how time is handled in the model, we created a starting date and tell the model how much time each 'step' should encompass 
 * further I create a variable counting the amount of weeks the model is running and a variable which returns the current week of the year
 */
	float step <- 0.5 #d; // this is the step, indicating the model is run every 0.5 day
 date starting_date <- date("2023-01-01");
	//date end_date <- date("2023-01-01");  //optional for a stop clause 
 int current_week <- current_date.week_of_year update: current_date.week_of_year;
	float weeks_passed;

	reflex weeks_busy {
		weeks_passed <- (current_date - starting_date) / #week; // Convert the difference (in seconds) to weeks

	}

	/* We continue with initilizing the model, we import all the files we will use in the model, and which the model will need to read
 * to create agents or grids
  */
	grid_file suitability_raster <- grid_file("../includes/suitability_utm.tif");
	file roads <- file("../includes/roads/motorways_links.shp");
	file wrong_study_area <- file("../includes/StudyAreaShp/Study_area_temp.shp");
	file boarlocs <- file("../includes/boar_locs_filtered.shp");
	file dens_kernel <- file("../includes/kernel5.tif");
	file faunapassage_file <- file("../includes/Faunapassages.shp");

	/* Setting the shape parameter, to indicate what the study area is, the world will be build ontop of this 'shape' geometry
	 */
	geometry shape <- envelope(suitability_raster);
	geometry road_graph <- envelope(roads);

	/* Initilizing input parameters for initial population size, initial age distribution and reproduction for both adult and yearling female boars  */
	int nb_sounder_init <- 10;
	float prob_repr_y <- 0.00989011; // Reproduction probability per time step
 float prob_repr_a <- 0.00989011; // Reproduction probability per time step
 list<float>
	age_distribution_prob_adult <- [0.39, 0.24, 0.15, 0.09, 0.06, 0.03, 0.02, 0.01, 0.01];

	/* before we start creating the model environment we will create TRUE or FALSE parameters for either good, bad or normal mast years
	 * At the start of each year it will be randomly decided (with a probability) if the year is good, bad or 'normal'
	 */
	bool bad_mast_year;
	bool good_mast_year;
	bool action_done_for_year <- false;

	reflex reset_action_flag when: (current_date - starting_date) mod 1 #y = 0 {
		action_done_for_year <- false;
		write ("action_done_for_year" + action_done_for_year);
	}

	reflex mast_year when: (current_week = 1 or current_week = 52) and not action_done_for_year {
		if rnd(10) < 2 and current_week = 1 {
			good_mast_year <- true;
			write ("good_mast" + good_mast_year);
		} else if rnd(10) < 2 and current_week = 1 {
			bad_mast_year <- true;
			write ("bad_mast" + bad_mast_year);
		} else {
			good_mast_year <- false;
			bad_mast_year <- false;
		}

		action_done_for_year <- true;
	}

	/* The following code blocks are initilzing the agents which exist in the model environment which the wild boars will interact with */
	grid kernel file: dens_kernel {
		rgb color <- hsb(grid_value / 10, 1.0, 1.0);
		geometry cell_shape <- square(1);
	}

	grid suitability_cell file: suitability_raster neighbors: 8 {
		bool contains_road <- false;
		bool contains_faunapassage <- false;
		rgb cell_color <- rgb(255, 255, 255); // This is white
 float cell_size <- 1000.0 #m;
		geometry cell_shape <- square(1);
		int boar_time <- 0;
		rgb color <- rgb(255, 255, 255); // default white color
 aspect default {
			draw shape color: #transparent border: #black;
		}

		/*I create a different aspect for the output which shows the utilization density, 
		 * rather than the suitability value */
		aspect density {
			draw shape color: boar_time_color border: #transparent;
		}

		rgb boar_time_color <- rgb(255, 255, 255);
		/* this is purely here for the display inside GAMA, the output of this is not used */
		reflex update_color {
			float intensity <- min(float(boar_time) / 100, 0.5);
			int color_shift <- round(255 * intensity);
			boar_time_color <- rgb(255 - color_shift, 255, 255 - color_shift);
		}

		init {
			color <- hsb(grid_value, 1.0, 5.0);
		}

		/* Important code which gives each cell a list of neighbors, 
		 * 2,3,4 respectively, which can then be called for when asking a specific cell
		 */
		list<suitability_cell> neighbors2 <- (self neighbors_at 2);
		list<suitability_cell> neighbors3 <- (self neighbors_at 3);
		list<suitability_cell> neighbors4 <- (self neighbors_at 4);
	}

	/* after creating the environmental agents we will now create those agents in the initialization of global
  * We will create all these agents ontop of the grid we created above
  * wildlife crossings
  * roads
  * initial boar locs (gbif sightings)
  * Sounders (+micro species)
  * adult wild boars
  */
	init {
		create faunapassages from: faunapassage_file {
			shape <- self.shape;
		}
		/* this code creates the roads, and also the inverse, which is not actually used in this model */
		create road from: roads {
			switch direction {
				match 0 {
					color <- #green;
				}

				match 1 {
					color <- #red;
					shape <- polyline(reverse(shape.points));
				}

				match 2 {
					color <- #blue;
					create road {
						shape <- polyline(reverse(myself.shape.points));
						direction <- 2;
						color <- #blue;
					}

				}

			}

		}

		/* after creating the road we want to create a variable which tells the suitability cell if it has a road inside or not
    * this parameter is initialized in suitability cell (bool contains_road) and is given a value true or false in this initilization code
    *The exact same happens for wildlife crossings below
    */
		loop a_cell over: suitability_cell {
			list<road> road_inside <- agents_overlapping(a_cell) of_species road;
			if length(road_inside) >= 1 {
				a_cell.contains_road <- true;
			}

		}

		loop a_cell over: suitability_cell {
			list<faunapassages> contains_faunapassage <- agents_overlapping(a_cell) of_species faunapassages;
			if length(contains_faunapassage) >= 1 {
				a_cell.contains_faunapassage <- true;
			}

		}

		create boar_locs from: boarlocs {
			shape <- self.shape;
		}

		/* Here we created the sounders, on random GBIF locations, as long as the suitbaility value is not nil (nil = nil does not return true) */
		create sounder number: nb_sounder_init {
			boar_locs chosen_boar_loc <- one_of(boar_locs);

			// Use the shape attribute of the chosen boar_loc for overlap
 my_cell <- one_of(suitability_cell where (each.shape overlaps chosen_boar_loc.shape and
			(each.grid_value = each.grid_value)));
			// Link the suitability cell of the sounder to a kernel cell
 loop each_kernel_cell over: kernel {
				if (my_cell.shape overlaps each_kernel_cell.shape) and !nil {
					my_kernel_cell <- each_kernel_cell;
					break; // Stop the loop once a matching kernel cell is found
 }

			}

			if (my_cell != nil) { // Just a check to make sure a cell was found
 begin <- my_cell;
				location <- my_cell.location;
			}
			/* In the following blocks we create micro species within the sounder to form a group. 
			 * The size of the group is determined by a max value and the density grid
			 * The initial age of each species is determined by putting the age_distribution_probability in the roulette wheel (assigned later)
			 */
			create adult_male number: round(7 * (my_kernel_cell.grid_value / 10)) {
				location <- myself.location; // set the sounder for this agent
 int index_m <- roulette_wheel(age_distribution_prob_adult);
				birth_date <- current_date - (2 #y + (rnd(index_m + 0.6, index_m + 0.9))) #y;
				my_cell <- myself.my_cell;
			}

			// Create initial adult females for this sounder:
 create adult_female number: round(8 * (my_kernel_cell.grid_value / 10)) {
				my_sounder <- myself;
				int index_f <- roulette_wheel(age_distribution_prob_adult);
				birth_date <- current_date - (2 #y + (rnd(index_f + 0.6, index_f + 0.9)) #y);
			}

			create male_yearling number: round(5 * (my_kernel_cell.grid_value / 10)) {
				my_sounder <- myself;
				birth_date <- current_date - (rnd(540, 690) #d);
			}

			create female_yearling number: round(8 * (my_kernel_cell.grid_value / 10)) {
				my_sounder <- myself;
				birth_date <- current_date - (rnd(540, 690) #d);
			}

			create piglet number: 30 * (my_kernel_cell.grid_value / 10) {
				my_sounder <- myself;
				birth_date <- current_date - (rnd(220, 330) #d);
			}

		}

	}

}

/* Now that we did all initializations we will describe the model, and the agents within the model - 
 * We start with creating a generic species, as seen below, which has attributes shared among all boars,
 * this includes things like age, birth_date and attributes for reflexes
 * Every reflex defined in this species can be used in all child species, and every attribute defined here 
 * will automatically also be defined in the child species unless specified differently
 * For example the default for the move comment = staying in its current cell, this is changed for the species that do move
 * Things like birth_date is unasigned, and a value is assigned in the specific species
 * A reflex like update_age, happens for each species that has an age assigned, a sounder for example does not have an age assigned, 
 * but the microspecies inside do
 */
species generic_species {
	list<float> probabilities;
	suitability_cell begin;
	suitability_cell my_cell;
	kernel my_kernel_cell;
	image_file my_icon;
	date birth_date;
	float age;

	reflex update_age when: every(1 #d) {
		if (birth_date != nil) {
			age <- (current_date - self.birth_date) / #year;
			// the /#year makes it so that the days (in which current_date and birth_date are calculated be divided so it ends up in years)
 }

	}

	/* Function for calculating probabilities 
	 * Returns a semi random index value between 0 and the 'length of the list you provide it with' depending on the values in said list
	 */
	int roulette_wheel (list<float> prob_list) {
		float r <- rnd(1.0); // random value between 0 and 1
 float cumulative_probability <- 0.0;
		loop i from: 0 to: (length(prob_list) - 1) {
			float prob <- prob_list[i];

			// Check for NaN and replace with 0
 if (prob != prob) {
				prob <- 0.0;
			}

			cumulative_probability <- cumulative_probability + prob;
			if (r <= cumulative_probability) {
				return i;
			}

		}

		return length(prob_list) - 1;
	}

	/* creating the attribute choose_cell, which represents a suitability_cell but is assigned nil on default */
	suitability_cell choose_cell {
		return nil;
	}
	/* Here make the current agent move towards the result of choose_cell, if choose_cell is not defined it looks if my_cell is nil or not,
 * if it is not nil the agent stays at the current position
 * additionally this is the place we keep track of boar_time, the unit for utilization density, every time a boar does move, it adds a boar_time to the location
 */
	reflex basic_move {
		my_cell <- choose_cell();
		if (my_cell != nil) {
			my_cell.boar_time <- my_cell.boar_time + 1;
			location <- my_cell.location;
		}

	}

}
/* Here we create the sounder agent, which consists of 4 micro species, additionally it has some attributes on its own, 
  * like the number of each microspecies inside it (local_*_*), size on the display (this is equal to its home range), a reflex for when it should die
  * choose_cell is chosing its current cell, making sure it does not move,
  * It can disperse its microspecies, or more precisely its yearlings, when they have not already dispersed this year
  * and lastly each sounder has a max group size defined as carrying capacity
  * There is also a code block for hunting, however, that is depricated for now. 
  */
species sounder parent: generic_species {
	int local_female_adults;
	int local_male_yearlings;
	int local_female_yearlings;
	int local_piglets;
	float size <- 4000.0;
	rgb color <- #black;

	aspect base {
		draw circle(size) color: color;
	}

	reflex update_local_counts {
		local_female_adults <- length(adult_female);
		local_male_yearlings <- length(male_yearling);
		local_female_yearlings <- length(female_yearling);
		local_piglets <- length(piglet);
	}

	reflex sounder_mortality {
		if length(adult_female + female_yearling) < 1 {
			do die;
		}

	}

	suitability_cell choose_cell {
		return my_cell;
	}

	species piglet parent: generic_species parallel: true {
		sounder my_sounder;
		/* quick tip - myself refers to the agent 2 tiers above, whereas self refers to the current agent 
 * in the following code self.age refers to the age of the piglet - we then create a female_yearling, everything inside that codeblock self would 
 * refer to the female yearling, thus we use myself, meaning the piglet once again. You can often hover over these to see what the model is refering to
 
 The rnd(10) <= 5 simply picks a random number between 0-10 and looks if it is smaller or equal to 5, if true the code gets executed, else the ELSE code gets executed.
 */
		reflex age {
			if self.age >= 1 {
				if rnd(10) <= 5 {
					create female_yearling number: 1 {
						my_sounder <- myself.my_sounder;
						birth_date <- myself.birth_date;
					}

					do die;
				} else {
					create male_yearling number: 1 {
						my_sounder <- myself.my_sounder;
						birth_date <- myself.birth_date;
					}

					do die;
				}

			}

		}

		/* here we code mortality, where we flip the daily death rate for each agent - if true the individual agent dies - good and bad mast year change the survival rate */
		reflex mortality when: every(1 #d) { // the when: statement makes sure this reflex is only done every day, instead of 0.5
 if flip(0.002328767) and not good_mast_year and not
			bad_mast_year {
				do die;
			} else if flip(0.002328767 * 2) and bad_mast_year {
				do die;
			} else if flip(0.002328767 / 2) {
				do die;
			} } }

	species female_yearling parent: generic_species parallel: true {
		sounder my_sounder;

		reflex mortality when: every(1 #d) {
			if flip(0.001726027) {
				do die;
			}

		}

		reflex age {
			if self.age >= 2 {
				create adult_female number: 1 {
					my_sounder <- myself.my_sounder;
					birth_date <- myself.birth_date;
				}

				do die;
			}

		}

		/* In the following code block we code reproduction - we do this by first creating an empty attribute representing the number of ofssprings,
 * we assign a pregnant interval simply to make sure they can only give birth once a year. 
 * 
 * The reproduce reflex is only called when the flip of the reproduction chance (daily chance to reproduce, making sure aorund 90% of the females reproduce per year)
 * and the time of the year is between march and june. 
 * The number of piglets created is dependend on a probability distribution - this probability distribution gets put into the roulette wheel to decide
 * how many offspring we get. Afterwards there is a chance this probability is increase if there is a good mast year. We assign the new piglet to the same sounder
 * as the female and assign the current birth_date to today. Further we make sure that the nb_offspring and preg interval are reset. 
 */
		int nb_offspring;
		int preg_interval <- 180 update: preg_interval + 1;

		reflex reproduce when: /* Create reproduction reflex */
		flip(prob_repr_a) and (current_week > 9) and (current_week < 22) and (self.preg_interval > 224) {
			list<float> prob_dist_sow_nb <- [0.01, 0.07, 0.16, 0.25, 0.25, 0.16, 0.07, 0.2, 0.01];
			int index_repro <- roulette_wheel(prob_dist_sow_nb);
			write ("nbpigs" + index_repro);
			if good_mast_year and flip(0.2) {
				nb_offspring <- round(index_repro * 1.5);
			} else {
				nb_offspring <- index_repro;
			}

			//if local_piglets > 50 {nb_offspring <- 0;} 
 create piglet number: nb_offspring {
				my_sounder <- myself.my_sounder;
				birth_date <- current_date;
			}

			set nb_offspring <- 0;
			set preg_interval <- 0;
		}

	}

	/* nothing new happens here */
	species male_yearling parent: generic_species parallel: true {
		sounder my_sounder;

		reflex age {
			if self.age >= 2 {
				create adult_male number: 1 {
					location <- myself.my_sounder.location;
					birth_date <- myself.birth_date;
					my_cell <- myself.my_sounder.my_cell;
				}

				do die;
			}

		}

		reflex mortality when: every(1 #d) {
			if flip(0.002219178) {
				do die;
			}

		}

	}

	/* The next code is for adult females, where we add something new
 * We make sure that when adults get older, their chance to die increase and their age cannot get above 14, the rest of the code is the same as female yearlings
 */
	species adult_female parent: generic_species parallel: true {
		int nb_offspring;
		sounder my_sounder;
		//int age <- round(current_date - birth_date) update: round(current_date - birth_date);
 reflex death_rate when: every(1 #d) {
			if flip(0.001) {
				do die;
			}

			if self.age > 10 and flip(0.001) { //double the chance to die if above the age of 10
 do die;
			}

			if self.age > 14 {
				do die;
			}

		}

		int preg_interval <- 180 update: preg_interval + 1;

		reflex reproduce when: /* Create reproduction reflex */
		flip(prob_repr_a) and (current_week > 9) and (current_week < 22) and (self.preg_interval > 224) {
			list<float> prob_dist_sow_nb <- [0.01, 0.07, 0.16, 0.25, 0.25, 0.16, 0.07, 0.2, 0.01];
			int index_repro <- roulette_wheel(prob_dist_sow_nb);
			if good_mast_year and flip(0.2) {
				nb_offspring <- round(index_repro * 1.5);
			} else if bad_mast_year and flip(0.2) {
				nb_offspring <- round(index_repro * 0.5);
			} else {
				nb_offspring <- index_repro;
			}

			//if local_piglets > 50 {nb_offspring <- 0;} 
 create piglet number: nb_offspring {
				my_sounder <- myself.my_sounder;
				birth_date <- current_date;
			}

			set nb_offspring <- 0;
			set preg_interval <- 0;
		}

	}

	/* in the following few code blocks we develope the dispersing aspect of the model
	 * Yearlings can disperse under certain conditions, which means they leave the sounder, and become a real agent and not just 
	 * a micro species. 
	 * Females are limited to have 2 yearlings disperse per year, whereas male yearlings can desperse in higher numbers. 
	 * *NOTE* we use myself.my_sounder.xx making - we are asking self.female_yearlings to create dispersing_females,
	 *  which then need to get parameters from the sounder of female_yearling. Myself thus refers to female_yearling, asking for information from their sounder
	 * Once the dispersing female is created, we ask the female yearling to die, so we don't duplicate our yearlings in the model. 
	 */
	float recently_female_dispersed <- 0.0 update: recently_female_dispersed + 0.5;

	reflex dispersing_female {
	// Ensure current cycle is week 28
 if ((current_week >= 25 and current_week <= 35) and (local_female_yearlings >= 4) and recently_female_dispersed < 50 #day) {
			if flip(0.2) {
				ask 2 among self.female_yearling {
					create dispersing_females number: 1 {
					//location <- self.location;  // 'self' refers to the sounder agent running this reflex.
 original_location <- myself.my_sounder.begin;
						my_cell <- myself.my_sounder.begin;
						birth_date <- myself.birth_date;
					}

					do die;
				}

			}

		}

	}

	float recently_male_dispersed <- 0.0 update: recently_male_dispersed + 0.5;

	reflex dispersing_male {
		if ((current_week >= 25 and current_week <= 30) and (local_male_yearlings >= 2) and recently_male_dispersed < 50 #day) {
			if flip(0.2) {
				ask (local_male_yearlings - 1) among self.male_yearling {
					create dispersing_males number: 1 {
						original_location <- myself.my_sounder.begin; // Assuming you have `my_cell` attribute in sounder to know its current cell.
 my_cell <- myself.my_sounder.begin;
						birth_date <- myself.birth_date;
					}

					do die;
				}

			}

		}

	}

	reflex Carrying_capacity {
	// Parameters
 int K <- 20; // Carrying capacity
 int N_0 <- length(adult_female + female_yearling + male_yearling); // Initial count
 int to_hunt <- N_0 - K;
		// Ensure no negative values
 if (to_hunt < 0) {
			to_hunt <- 0;
		}

		// Randomly select and "hunt" the calculated number of yearlings/adults
 if flip(0.1 * to_hunt) {
			ask rnd(1, to_hunt) among (adult_female + female_yearling + male_yearling) {
				do die;
			}

		}

	}

	/* bool is_hunting_activated <- false;

	reflex hunting_activation when: nb_sounders > 80 and not is_hunting_activated {
		is_hunting_activated <- true;
		write ("hunting" + is_hunting_activated);
	}

	reflex hunting_sounders when: is_hunting_activated {
		if (current_week >= 38) {
			if flip(0.005) {
				ask 1 among sounder {
					do die;
				}

			}

		} else if (current_week > 51) or nb_sounders < 30 {
		// Reset the activation for the next year
			is_hunting_activated <- false;
		}

	} */
}

species adult_male parent: generic_species skills: [moving] parallel: true {
	float size <- 3000.0;
	rgb color <- #yellow;

	aspect base {
		draw circle(size) color: color;
	}

	suitability_cell choose_cell {
		if flip(0.5) {
			return my_cell;
		} else {
			probabilities <- [];
			list<suitability_cell> neighbors <- [my_cell];
			list<suitability_cell> potential_neighbors <- my_cell.neighbors;
			loop cell over: potential_neighbors {
				if (cell.grid_value != cell.grid_value) { // Check for NaN and direct road/fence
 continue; // Skip this cell
 }

				if cell.contains_road = true and not cell.contains_faunapassage {
					if flip(0.8) {
						if flip(0.5) {
							do die;
						} else {
							continue;
						}

					}

				}

				// Check the path between my_cell and the potential destination
 neighbors <- neighbors + cell;
			}

			// Construct probabilities list
 float total_suitability <- 0.0;
			loop cell over: neighbors {
				total_suitability <- total_suitability + cell.grid_value;
			}

			loop cell over: neighbors {
				float probability <- cell.grid_value / total_suitability;
				probabilities <- probabilities + [probability];
			}

			int index <- roulette_wheel(probabilities);
			return neighbors[index];
		}

	}

	reflex mortality when: every(1 #d) {
		if flip(0.001369863) {
			do die;
		}

		if self.age > 8 and flip(0.001369863) {
			do die;
		}

		if self.age > 14 {
			do die;
		}

	}

}

species dispersing_males parent: generic_species skills: [moving] parallel: true {
	float size <- 3000.0;
	rgb color <- #purple;

	aspect base {
		draw circle(size) color: color;
	}

	suitability_cell original_location;
	suitability_cell choose_cell {
		if flip(0.5) {
			return my_cell;
		} else {
			probabilities <- [];
			list<suitability_cell> neighbors <- [my_cell];
			list<suitability_cell> potential_neighbors <- my_cell.neighbors;
			loop cell over: potential_neighbors {
				if (cell.grid_value != cell.grid_value) { // Check for NaN and direct road/fence
 continue; // Skip this cell
 }

				if cell.contains_road = true and not cell.contains_faunapassage {
					if flip(0.8) {
						continue;
					}

				}

				// Check the path between my_cell and the potential destination
 neighbors <- neighbors + cell;
			}

			// Construct probabilities list
 float total_suitability <- 0.0;
			loop cell over: neighbors {
				total_suitability <- total_suitability + cell.grid_value;
			}

			loop cell over: neighbors {
				float probability <- cell.grid_value / total_suitability;
				probabilities <- probabilities + [probability];
			}

			int index <- roulette_wheel(probabilities);
			return neighbors[index];
		}

	}

	reflex age {
		if self.age >= 2 {
			create adult_male {
				location <- myself.location;
				birth_date <- myself.birth_date;
				my_cell <- myself.my_cell;
			}

			do die;
		}

	}

	reflex do_die when: every(1 #d) {
		if flip(0.002219178) {
			do die;
		}

	}

}

species dispersing_females parent: generic_species skills: [moving] {
	float size <- 3000.0;
	rgb color <- #red;

	aspect base {
		draw circle(size) color: color;
	}

	suitability_cell original_location;
	int weeks_moved <- 0;

	reflex new_sounder {
		bool no_sounders_nearby <- true;
		if !nil and my_cell.grid_value > 0.23 {
			loop each_cell over: my_cell.neighbors3 {
				if (length(sounder inside each_cell) > 0) {
					no_sounders_nearby <- false;
					break; // exit the loop if a sounder is found in any neighboring cell
 }

			}

			if (no_sounders_nearby) and flip(0.1) {
				create sounder {
					birth_date <- myself.birth_date;
					location <- myself.location;
					my_cell <- myself.my_cell;
					begin <- my_cell;
					create female_yearling number: 3 {
						my_sounder <- myself;
						birth_date <- myself.birth_date;
					}

				}

				do die;
				ask 2 among (self.dispersing_female where distance(self.location, myself.location) <= 15 * 1000) {
					do die;
				}

			}

		}

	}

	suitability_cell choose_cell {
		if flip(0.5) {
			return my_cell;
		} else {
			probabilities <- [];
			list<suitability_cell> neighbors <- [my_cell];
			list<suitability_cell> potential_neighbors <- my_cell.neighbors;
			list<suitability_cell> valid_neighbors <- [];

			// write("neighbors" +my_cell.neighbors2);
 // Filter out any neighbors with NaN grid values
 loop cell over: potential_neighbors {
			// Check if cell is within a maximum of 10 cells away from original location
 if (bool((cell neighbors_at 2) where
				(max(abs(cell.grid_x - original_location.grid_x), abs(cell.grid_y - original_location.grid_y)) <= 10))) {
					valid_neighbors <- valid_neighbors + cell;
				}

			}
			// If no valid neighbors, return nil
 if (empty(valid_neighbors)) {
				return my_cell;
			} else {
				loop cell over: valid_neighbors {
					if (cell.grid_value != cell.grid_value) { // Check for NaN and direct road/fence
 continue; // Skip this cell
 }

					if cell.contains_road = true and not cell.contains_faunapassage {
						if flip(0.8) {
							continue;
						}

					}

					// Check the path between my_cell and the potential destination
 neighbors <- neighbors + cell;
				}

				// Construct probabilities list
 float total_suitability <- 0.0;
				loop cell over: neighbors {
					total_suitability <- total_suitability + cell.grid_value;
				}

				loop cell over: neighbors {
					float probability <- cell.grid_value / total_suitability;
					probabilities <- probabilities + [probability];
				}

				int index <- roulette_wheel(probabilities);
				return neighbors[index];
			}

		}

	}

	reflex do_die when: every(1 #d) {
		if flip(0.001205479) {
			do die;
		}

	}

}

species road {
	int direction;
	rgb color;

	aspect geom {
		draw shape color: color;
	}

}

species faunapassages {

	aspect geom {
		draw shape color: #green;
	}

}

species boar_locs {

	aspect geom {
		draw shape color: #red;
	}

}

experiment sounder_movement type: gui {
	parameter "Initial sounders:" var: nb_sounder_init min: 1 max: 1000 category: sounder;
	output {
		display main_display type: java2D {
			grid suitability_cell border: #gray;
			species road aspect: geom;
			species sounder aspect: base;
			species dispersing_females aspect: base;
			species dispersing_males aspect: base;
			species adult_male aspect: base;
			//species boar_locs aspect: geom;	
 }

		display density {
			species suitability_cell aspect: density;
			species road aspect: geom;
			//species sounder aspect: base;
 //species dispersing_females aspect: base;
 //species dispersing_males aspect: base;
 //species adult_male aspect: base;

		}

		display Population_dynamics refresh: every(2 #day) {
			chart "Population dynamics overview: time series" type: series x_tick_unit: 180 x_serie_labels: string(current_date, "MM/yy") background: rgb(47, 47, 47) color: #white
			position: {0, 0} size: {1, 1} x_label: "Time (weeks)" y_label: "Abundance (nr)" {
				data "total population sounders" value: nb_sounders color: #pink thickness: 3 marker: false;
				data "total population size" value: (nb_male_adults + nb_female_adults + nb_male_yearlings + nb_female_yearlings + nb_piglets) color: #white marker: false thickness: 3;
				data "adults (nr)" value: (nb_female_adults + nb_male_adults) color: rgb(46, 204, 113) marker: false thickness: 2;
				data "yearlings (nr)" value: (nb_female_yearlings + nb_male_yearlings + nb_dispersing_females + nb_dispersing_males) color: rgb(231, 76, 60) marker: false thickness: 2;
				data "piglets (nr)" value: nb_piglets color: rgb(52, 152, 219) marker: false thickness: 2;
			}

		}

	}

}

