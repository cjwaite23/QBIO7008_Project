# QBIO7008_Project

# Author: Callum Waite

Research Supervisor:

-   Dr Simon Hart, School of Biological Sciences, University of Queensland

**Project and Repository Description**

This is a repository containing most of the code and figures I produced and used during my QBIO7008 Masters research placement with Simon Hart. The purpose of this project was to explore and fit state-space analyses to a large capture-recapture dataset of Trinidadian guppies (*Poecilia reticulata*).

Due to the size of the raw data (obtained from Griffiths *et al*., (2020)), it is not present in that form here in this repository. The same applies for any large data files that were generated in one of the scripts but are too great in size to be stored on GitHub. For this reason, some file-save or -load locations will reference my personal hard-drive. This applies especially to large manipulated data files and Bayesian analyses that are created with many iterations. All generated figures in the code are stored as pdf vector files.

**General Code Structure**

In general, R scripts in this repository are structured similarly. They all begin with a description that outlines what the general function of that script is. Often the scripts will be quite long, which is why they are followed by a loading of libraries section, and then potentially many more sections, all labelled to describe their exact function. In some instances, similar code is provided many times in the one script with slight changes made to accomodate small changes in data requirements. Not all of these repeat instances will be adequately commented.

Important notes:

-   the structure of most folders (Analyses, Results, Figures) are in the form of two sub-folders - *MARK_Analyses* and *SSMs*, to represent the two components of analysis performed here. Within *SSMs* we have another split into *jolly-seber* and *pop_growth* to denote files related to different types of state-space models that I fit

-   The state-space models fitted under a Jolly-Seber capture-recapture formulation were too computationally intensive to be run to convergence. No results (except one 100 iteration run) and no figures are therefore provided.

-   This repository does not contain a thorough catalogue of every possible idea I came up with throughout this project, because I generated a lot of rough working scripts that are often repeated to near exactness by scripts that do lie in this repository. Please do not hesitate to contact me if you are interested in further ideas and conclusions or have queries about any scripts.

**Repository Contents**

***1_Data***

-   No raw data was actually small enough in size to be included in this repository. The only raw data is a file directly obtained from Griffiths *et al.* (2020), and climatic data obtained from WorldClim.org. Both these files are located on my personal hard-drive and can be provided upon request, however they are both publicly available.

-   The raw capture-recapture data file is a .csv file named after the paper by Griffiths *et al.* (2020), *"Individual_differences_determine_the_strength_of_ecological_interactions.csv"*. It contains detailed C-R presence-absence data for every individual guppy sampled in all four streams throughout the 97 sampling occasions. It also contains weight and size measurements, of which I did not analyse.

-   The climatic data contains global precipitation, max temp, and min temp data at 2.5' precision. They are in the form of .tiff forms for each month from 2000-2018.

***2_Data_manipulation***

-   "*MARK_encounter_histories_format.R*": R script to take raw C-R data from my personal hard-drive and convert mark-recapture data of guppies into .inp files for analysis with program MARK and into R data frames for analysis with SSMs

-   *"SSM_capture_histories.R"*: R script to generates capture histories for each river to be analysed by state-space formulations.

-   *"worldclim_data_extraction.R*": R script to extract relevant data from the raw WorldClim data. This means filtering it to the relevant location of the sites in Northern Trinidad, and to the relevant months (2008-2016).

-   FOLDER "*data_files"*: manipulated data files

    -   FOLDER *"MARK_encounter_histories"*: folder containing all manipulated files for each stream in the form ready to be analysed by passing to MARK

    -   *"rmark_histories.RData"*: data-frames produced by *"MARK_encounter_histories_format.R"*

    -   *"capture-histories_df.RData"*: capture history dataframes for each sex and stream produced by *"SSM_capture_histories.R"*

    -   *"climate_data.RData"*: manipulated climatic data to only specify for the stream locations and the months in which sampling occurred.

***3_Analyses***

-   FOLDER *"MARK_analyses"*: folder containing analyses pertaining to those models implemented in program MARK.

    -   *"rmark_best_js_models.R"*: R script to implement and select the best Jolly-Seber models in MARK by interfacing through package "rmark".

-   FOLDER *"SSMs"*: folder containing analyses pertaining to the implementation of any state-space models

    -   FOLDER *"pop_growth"*: folder containing state-space analyses of pop size estimates generated by MARK. This is the folder containing the analyses that ultimately led to figures and tangible results

        -   FOLDERS "*jags_models"* and *"stan_models*": folders containing raw .jags or .stan files with the model descriptions for implementation in R. Specific scripts are named relating to their formulation - no subscript represents the simplest implementation, "*\_dd"* represents a density-dependent formulation

        -   "*SSM_pop_growth_final.R"*: largest R script, containing the full process of the successful state-space implementation in JAGS and linear regressive analysis of growth-rate. This is the final process that was followed to generate the results and figures shown in the report.

        -   *"SSM_pop_growth_stan_R*": R script to implement the Stan pop-growth state-space models

        -   *"kalman_filter_exploration.R"*: R script exploring the use of the Kalman filter from the package "dlm" for implementation on the data. Not reported on in the manuscript but does provide a good comparison to the final model used.

    -   FOLDER *"jolly_seber"*: folder containing **attempted** implementations of the direct state-space formulation of the Jolly-Seber capture-recapture model.

        -   FOLDERS "*jags_models"* and *"stan_models*": folders containing raw .jags or .stan files with the model descriptions for implementation in R. Specific scripts are named relating to their formulation - no subscript represents the simplest implementation, "*\_dd"* represents a density-dependent formulation

        -   FOLDER *"hpc_code*": contains scripts and bash-scripts for implementing these JAGS models within UQ's QBIO HPC cluster. Again, files are named by their formulation, with a subscript of *"b0"* indicating the use of just a constant survival parameter, while *"b0b1N*" indicates the inclusion of survival being density dependent with coefficient b1. BAsh scripts have a suffix *"bs\_*"

        -   *"js_ssm\_(.)\_(dd).R*": R scripts to implement the jags or stan formulations of the Jolly-Seber state-space modelm with (dd) or without density-dependence.

        -   *"moving_window_approach.R*": R script to investigate the use of the moving window approach I developed and conceptualised in the manuscript. Not successfully run in JAGS due to language constraints, but the premise remains relevant

***4_Results***

-   FOLDER *"MARK_analyses"*: folder containing results pertaining to those models implemented in program MARK.

    -   FOLDER *"rmark_best_models"*: folder containing the results from the best models generated by rmark

        -   *"(stream suffix)(sex)\_best_model.RData":* files containing the best MARK Jolly-Seber models (always time-dependent for all parameters) for each stream and sex combination.

        -   *"(stream)\_B*\_*derived.RData"*: files containing dataframes of the best models and their derived estimates of parameters including survival, recruitment and capture probabilities, population sizes and number of recruits. Each variable has a mean and a 95% confidence range.

-   FOLDER *"SSMs"*: folder containing results pertaining to the implementation of any state-space models

    -   FOLDER *"pop_growth"*: folder containing state-space results of pop size estimates generated by MARK. This is the folder containing the key figures and tangible results

        -   FOLDER *"final_riv_sex_models"*: folder containing the results relating to the final river and sex models

            -   *"ssm_preds.RData"*: file containing the predicted population sizes and all other estimated variables for each time period and each sex in each stream, as generated by the state space analyses of MARK results

            -   *"reg_data.RData"*: file containing the transformed dataframe of results prepared for regressive analysis.

    -   FOLDER *"jolly_seber"*: folder containing results of a single attempted implementation of the direct state-space formulation of the Jolly-Seber capture-recapture model.

        -   *"CA_js_ssm_phit.RData*": example data of an implementation of the state-space formulation of the Jolly-Seber model. The data run was the full CA (Caigual) dataset, with survival (phi) only dependent on time (t). Only 100 iterations were run on the HPC and this still took 178 minutes.

***5_Figures***

-   FOLDER *"SSMs"*: folder containing the figures produced by the state_space analyses of the data
    -   FOLDER "*pop_growth*": only pop_growth SSMs had any interpretable results

        -   "*(all/CA/LL/TY/UL)\_plot.pdf/jpg"*: figures from SSM analysis comparing the population estimates at each stream (or all of them) against the observed values.

        -   *"visreg_plot.pdf/jpg"*: partial residual plots of key predictors against growth rate for the final regressive model produced.

Note that no manuscripts or slides are provided in this repository.

Have a nice day and thanks for stopping by!
