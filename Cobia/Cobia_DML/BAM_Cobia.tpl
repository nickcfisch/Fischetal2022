//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##
//##  SEDAR 58  SA Cobia, 2019
//##  (Modified from: SEDAR 50  SA BLT assessment model, 2017)
//##  NMFS, Beaufort Lab, Sustainable Fisheries Branch
//##
//## #Version with Thorson Dirichlet-Multinomial Linear 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
DATA_SECTION

!!cout << "Starting Beaufort Assessment Model" << endl;
!!cout << endl;
!!cout << "                BAM!" << endl;
!!cout << endl;

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: set-up section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

// Starting and ending years of the model (year data starts)
init_int styr;
init_int endyr;

// Starting and ending years to estimate recruitment deviations from S-R curve
init_int styr_rec_dev;
init_int endyr_rec_dev;
// Ending years of 3 phases of constraints on recruitment deviations
// (allows possible heavier constraint (weights defined later) in early and late period, with lighter constraint in the middle)
init_int endyr_rec_phase1;
init_int endyr_rec_phase2;
// ending years of selectivity block 1
//init_int endyr_selex_phase1_comm; // comm
init_int endyr_selex_phase1_GR; // GR
// ending years of selectivity block 2
//init_int endyr_selex_phase2_comm; // comm
//init_int endyr_selex_phase2_cL; // cL
//init_int endyr_selex_phase2_GR; // GR

//number assessment years
number nyrs;
number nyrs_rec;

//this section MUST BE INDENTED!!!
 LOCAL_CALCS
   nyrs=endyr-styr+1.;
   nyrs_rec=endyr_rec_dev-styr_rec_dev+1.;
 END_CALCS
 
// Number of ages in population model(classes are 1,...,N+; assumes last age is plus group
init_int nages;
// Vector of ages for age bins in population model, last is a plus group
init_vector agebins(1,nages);

//Total number of ages used to match age comps: plus group may differ from popn, first age must not
init_int nages_agec;

//Vector of ages for age bins in age comps
init_vector agebins_agec(1,nages_agec);

// Number length bins used to match length comps and width of bins (mm)
init_int nlenbins;          //used to match data
init_number lenbins_width;  //width of length bins (mm)

// Vector of length bins (mm; midpoint of bin) used to match length comps and bins used to compute plus group 
init_vector lenbins(1,nlenbins);
 
// Max value of F used in spr and msy calculations
init_number max_F_spr_msy;
// Number of iterations in spr calculations
init_int n_iter_spr;
//Total number of iterations for msy calcs
int n_iter_msy;
 LOCAL_CALCS
		n_iter_msy=n_iter_spr; 
 END_CALCS

// Starting and ending years to compute arithmetic average recruitment for SPR-related values
init_int styr_rec_spr;
init_int endyr_rec_spr;
//Arithmetic average recruitment for SPR-related values
number nyrs_rec_spr;
 LOCAL_CALCS
   nyrs_rec_spr=endyr_rec_spr-styr_rec_spr+1.;
 END_CALCS 

 
// Number of years at end of time series over which to average sector Fs, for weighted selectivities
init_int selpar_n_yrs_wgted;
// Multiplicative bias correction of recruitment (may set to 1.0 for none or negative to compute from recruitment variance)
init_number set_BiasCor;

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: observed data section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//######################################################################################################
//## Commercial Handline
//######################################################################################################
//## comm ######## INDEX ########
//## comm ## Starting and ending years of CPUE index
//init_int styr_comm_cpue;                                             
//init_int endyr_comm_cpue;
//## comm ## Observed index CPUE and CVs                                            
//init_vector obs_comm_cpue(styr_comm_cpue,endyr_comm_cpue);
//init_vector comm_cpue_cv(styr_comm_cpue,endyr_comm_cpue);

//## comm ######## LANDINGS ########
//## comm ## Starting and ending years of landings (bt50 landings+discards; includes comm and cO (commercial other))
init_int styr_comm_L;
init_int endyr_comm_L;
//## comm ## Observed landings (1000 lbs) and assumed CVs
init_vector obs_comm_L(styr_comm_L,endyr_comm_L);
init_vector comm_L_cv(styr_comm_L,endyr_comm_L);

//## comm ######## LENGTH COMPS ########
//## comm ## Number and vector of years of length compositions to be pooled
init_int nyr_comm_lenc_pool;
init_ivector yrs_comm_lenc_pool(1,nyr_comm_lenc_pool);
//# Annual sample size (nfish) of length comp data; used to weight years for pooling
init_vector nfish_comm_lenc_pool(1,nyr_comm_lenc_pool);
//## comm ## Number and vector of years of length compositions, after pooling
init_int nyr_comm_lenc;
init_ivector yrs_comm_lenc(1,nyr_comm_lenc);
//## comm ## Sample size of length comp data (first row observed n.trips, second row n.fish)
init_vector nsamp_comm_lenc(1,nyr_comm_lenc);
init_vector nfish_comm_lenc(1,nyr_comm_lenc);
//## comm ## Observed length comps (3cm bins; proportions by year)
init_matrix obs_comm_lenc(1,nyr_comm_lenc,1,nlenbins);

//######################################################################################################
//## Recreational Headboat
//######################################################################################################
//## HB ######## INDEX ########
//## HB ## Starting and ending years of CPUE index
init_int styr_HB_cpue;                                             
init_int endyr_HB_cpue;
//## HB ## Observed index CPUE and CVs                                            
init_vector obs_HB_cpue(styr_HB_cpue,endyr_HB_cpue);//Observed CPUE
init_vector HB_cpue_cv(styr_HB_cpue,endyr_HB_cpue); //CV of cpue

//######################################################################################################
//## General Recreational
//######################################################################################################
//## GR ######## LANDINGS ########
//## GR ## Starting and ending years of landings (bt50 landings+discards)
init_int styr_GR_L;
init_int endyr_GR_L;
//## GR ## Observed landings (1000 lbs) and assumed CVs
init_vector obs_GR_L(styr_GR_L,endyr_GR_L);   //vector of observed landings by year 
init_vector GR_L_cv(styr_GR_L,endyr_GR_L);    //vector of CV of landings by year

//## GR ######## LENGTH COMPS ########
//## GR ## Number and vector of years of length compositions
//init_int nyr_GR_lenc;
//init_ivector yrs_GR_lenc(1,nyr_GR_lenc);
//## GR ## Sample size of length comp data (first row observed n.trips, second row n.fish)
//nit_vector nsamp_GR_lenc(1,nyr_GR_lenc);
//init_vector nfish_GR_lenc(1,nyr_GR_lenc);
//## GR ## Observed length comps (3cm bins; proportions by year)
//init_matrix obs_GR_lenc(1,nyr_GR_lenc,1,nlenbins);

//  Age compositions
init_int nyr_GR_agec;
init_ivector yrs_GR_agec(1,nyr_GR_agec);
init_vector nsamp_GR_agec(1,nyr_GR_agec);
init_vector nfish_GR_agec(1,nyr_GR_agec);
init_matrix obs_GR_agec(1,nyr_GR_agec,1,nages_agec);

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: parameter section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//######################################################################################################
//## Parameter values and initial guesses
//######################################################################################################
//######## Population ########
init_vector set_Linf(1,7);				// VonBert Linf (mmFL)
init_vector set_K(1,7);					// VonBert K
init_vector set_t0(1,7);				// VonBert t0
init_vector set_len_cv(1,7);			// CV of length at age
//########Landings growth curve
init_vector set_Linf_L(1,7);            // VonBert Linf (mmFL)
init_vector set_K_L(1,7);               // VonBert K
init_vector set_t0_L(1,7);              // VonBert t0
init_vector set_len_cv_L(1,7);          // CV of length at age
//########Female only growth curve (popl and landings)
init_vector set_Linf_F(1,7);            // VonBert Linf (mmFL)
init_vector set_K_F(1,7);               // VonBert K
init_vector set_t0_F(1,7);              // VonBert t0
init_vector set_len_cv_F(1,7);          // CV of length at age
//######## Constant M ########
init_vector set_M_constant(1,7);		// constant M (used only to compute MSST=(1-M)SSBmsy)     
//######## StockRecruitment ########
init_vector set_steep(1,7);         	// SR steepness parameter
init_vector set_log_R0(1,7);        	// SR log_R0 parameter
init_vector set_R_autocorr(1,7);    	// SR recruitment autocorrelation (lag 1)
init_vector set_rec_sigma(1,7);     	// standard deviation of recruitment in log space
//######## DirichletMultinomial ########
init_vector set_log_dm_comm_lc(1,7);  	// Dirichlet-multinomial overdispersion parameter (log-space): comm length comps
//init_vector set_log_dm_cL_lc(1,7); 		// Dirichlet-multinomial overdispersion parameter (log-space): cL length comps
//init_vector set_log_dm_HB_lc(1,7);  	// Dirichlet-multinomial overdispersion parameter (log-space): HB length comps
//init_vector set_log_dm_GR_lc(1,7);  	// Dirichlet-multinomial overdispersion parameter (log-space): GR length comps
init_vector set_log_dm_GR_ac(1,7);    //Dirichlet-multinomial overdispersion parameter
//######## Selectivity ########
init_vector set_selpar_A50_comm1(1,7);	// comm age at 0.5 selectivity
init_vector set_selpar_slope_comm1(1,7);	// comm slope of ascending limb
//init_vector set_selpar_A50_comm2(1,7);	// comm age at 0.5 selectivity (block 2)
//init_vector set_selpar_slope_comm2(1,7);	// comm slope of ascending limb (block 2)
// init_vector set_selpar_A502_comm2(1,7);	// comm L502 (block 2)
// init_vector set_selpar_slope2_comm2(1,7);	// comm slope of descending limb (block 2)
//init_vector set_selpar_A50_comm3(1,7);	// comm age at 0.5 selectivity (block 3)
//init_vector set_selpar_slope_comm3(1,7);	// comm slope of ascending limb (block 3)
//init_vector set_selpar_A50_cL1(1,7);	// cL age at 0.5 selectivity
//init_vector set_selpar_slope_cL1(1,7);	// cL slope of ascending limb
//init_vector set_selpar_A50_cL2(1,7);	// cL age at 0.5 selectivity (block 2)
//init_vector set_selpar_slope_cL2(1,7);	// cL slope of ascending limb (block 2)
//// init_vector set_selpar_A502_cL2(1,7);	// cL L502 (block 2)
//// init_vector set_selpar_slope2_cL2(1,7);	// cL slope of descending limb (block 2)
//init_vector set_selpar_A50_cL3(1,7);	// cL age at 0.5 selectivity (block 3)
//init_vector set_selpar_slope_cL3(1,7);	// cL slope of ascending limb (block 3)
init_vector set_selpar_A50_GR1(1,7); 	// GR age at 0.5 selectivity
init_vector set_selpar_slope_GR1(1,7);	// GR slope of ascending limb
init_vector set_selpar_A50_GR2(1,7);	// GR age at 0.5 selectivity (block 2)
init_vector set_selpar_slope_GR2(1,7);	// GR slope of ascending limb (block 2)
// init_vector set_selpar_A502_GR2(1,7);	// GR L502 (block 2)
// init_vector set_selpar_slope2_GR2(1,7);	// GR slope of descending limb (block 2)
//init_vector set_selpar_A50_GR3(1,7);	// GR age at 0.5 selectivity (block 3)
//init_vector set_selpar_slope_GR3(1,7);	// GR slope of ascending limb (block 3)
//######## IndexCatchability ########
//init_vector set_log_q_comm(1,7);      	// comm CPUE (log q)
//init_vector set_log_q_cL(1,7);      	// cL CPUE (log q)
init_vector set_log_q_HB(1,7);      	// HB CPUE (log q)
//######## FishingMortality ########
init_vector set_F_init(1,7);  			// initial F (not log space) 
init_vector set_log_avg_F_comm(1,7);		// comm log mean F
//init_vector set_log_avg_F_cL(1,7);		// cL log mean F
init_vector set_log_avg_F_GR(1,7);		// GR log mean F

//######################################################################################################
//## Dev vectors
//######################################################################################################
init_vector set_log_F_dev_comm(1,3); 		// comm F devs
//init_vector set_log_F_dev_cL(1,3);		// cL F devs
init_vector set_log_F_dev_GR(1,3);		// GR F devs
init_vector set_log_RWq_dev(1,3);		// Random walk on q
init_vector set_log_rec_dev(1,3);		// recruitment devs
init_vector set_log_Nage_dev(1,3);		// Nage devs

//######## F dev initial guesses ########
//## comm (1962 - 2015)
init_vector set_log_F_dev_comm_vals(styr_comm_L,endyr_comm_L);
//## cL (1958 - 2015)
//init_vector set_log_F_dev_cL_vals(styr_cL_L,endyr_cL_L);
// ## GR (1973 - 2015)
init_vector set_log_F_dev_GR_vals(styr_GR_L,endyr_GR_L);

//######## Rec dev initial guesses (1958 - 2015) ########
init_vector set_log_rec_dev_vals(styr_rec_dev,endyr_rec_dev);

//######## initial N age devs, all ages but the first one (2 to 15) ########
init_vector set_log_Nage_dev_vals(2,nages);             

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: likelihood weights section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
init_number set_w_L;			// landings

//init_number set_w_I_comm;         // comm index
//init_number set_w_I_cL;         // cL index
init_number set_w_I_HB;         // HB index
init_number set_w_lc_comm;        // comm length comps
//init_number set_w_lc_cL;        // cL length comps
//init_number set_w_lc_GR;		// GR length comps
init_number set_w_ac_GR;		//weight for the Recreational age comps  
init_number set_w_Nage_init;    // log N.age.dev residuals (initial abundance)for fitting initial abundance at age (excluding first age)
init_number set_w_rec;          // SR residuals (for fitting SR curve)
init_number set_w_rec_early;    // constraint on early recruitment deviations
init_number set_w_rec_end;      // constraint on ending recruitment deviations
init_number set_w_fullF;        // penalty if F exceeds 3.0 (reduced by factor of 10 each phase, not applied in final phase of optimization) full F summed over fisheries
init_number set_w_Ftune;        // weight on tuning F (penalty not applied in final phase of optimization)

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: miscellaneous stuff section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

// Length-weight parameter a (W=aL^b); mm to kg
init_number wgtpar_a;
// Length-weight parameter b (W=aL^b); mm to kg
init_number wgtpar_b;

// Length-batchFecundity parameter a (BF=a + bL); mm to eggs
//init_number batchfecpar_a;
// Length-batchFecundity b (BF=a + bL); mm to eggs
//init_number batchfecpar_b;
// Number of batches spawned per mature female per year
//init_number nbatch;
// Value by which to scale fecundity unit (e.g. 1000000 to diplay fecundity in millions)
//init_number fecpar_scale;

// vector of maturity-at-age for females (ages 1 - 15 )
init_vector maturity_f_obs(1,nages);
// Proportion female by age (assumed 50:50 sex ratio)
init_vector prop_f_obs(1,nages);
// time of year (as fraction) for spawning
init_number spawn_time_frac;
// age-dependent natural mortality at age (ages 1 - 15 )
init_vector set_M(1,nages);

// Spawner-recruit parameters:  SR function switch (integer 1=Beverton-Holt, 2=Ricker)
init_int SR_switch;

// switch for rate increase in q: Integer value (choose estimation phase, negative value turns it off)
init_int set_q_rate_phase;
// annual positive rate of increase on all fishery dependent q's due to technology creep
init_number set_q_rate;

// density dependence on fishery catchability coefficients (DDq) switch: Integer value (choose estimation phase of random walk, negative value turns it off)
init_int set_q_DD_phase;
// q_DD exponent, value of zero is density independent, est range is (0.1,0.9)
init_number set_q_DD_beta;
// SE of q_DD exponent (0.128 provides 95% CI in range 0.5)
init_number set_q_DD_beta_se;
// Age to begin counting q_DD (should be age near full exploitation)
init_int set_q_DD_stage;      //age to begin counting biomass, should be near full exploitation

// Variance (sd^2) of fishery dependent random walk catchabilities (0.03 is near the sd=0.17 of Wilberg and Bence) 
init_number set_RWq_var;     //assumed variance of RW q

// Tuning F (not applied in last phase of optimization, or not applied at all if penalty weight=0)
init_number set_Ftune;
// Year for tuning F
init_int set_Ftune_yr;

// threshold sample sizes ntrips (>=)for length comps (set to 99999.0 if sel is fixed): 
init_number minSS_comm_lenc;	// comm len comps
//init_number minSS_cL_lenc;	// cL len comps
// init_number minSS_HB_lenc;	// HB len comps	(to be removed for bt50)
//init_number minSS_GR_lenc;	// GR len comps
//threshold sample sizes for age comps
init_number minSS_GR_agec;

// Input for deterministic F-based projections
// Last year of projections, must be later than assessment endyr by default
init_int endyr_proj;	// Projection end year (must be later than assessment endyr)
init_int styr_regs;  	// Apply current F until styr_regs, then the projection F
init_int Fproj_switch;  // Switching indicating value to use for defining projection F: 1=Fcurrent, 2=Fmsy, 3=F30, 4=F40
init_number Fproj_mult; // Multiplier 'c' applied to compute projection F, for example Fproj=cFmsy
// Calculate projection start year
int styr_proj;			
 LOCAL_CALCS
   styr_proj=endyr+1;
 END_CALCS
 
// Aging error matrix (columns are true age 1- 15 , rows are ages as read for age comps: columns should sum to one)
init_matrix age_error(1,nages,1,nages);

//------------------------------------<< 999 >>-------------------------------------
// END OF READING IN VALUES FROM .dat file
//------------------------------------<< 999 >>-------------------------------------


// #######Indexing integers for year(iyear), age(iage),length(ilen) ###############
int iyear;
int iage;
int ilen;
int ff;

number sqrt2pi;
number g2mt;                    //conversion of grams to metric tons 
number g2kg;                    //conversion of grams to kg   
number g2klb;                   //conversion of grams to 1000 lb   
number mt2klb;                  //conversion of metric tons to 1000 lb
number mt2lb;                   //conversion of metric tons to lb
number dzero;                   //small additive constant to prevent division by zero
number huge_number;             //huge number, to avoid irregular parameter space

init_number end_of_data_file;
//this section MUST BE INDENTED!!!
 LOCAL_CALCS
   if(end_of_data_file!=999)
   {
       cout << "*** WARNING: Data File NOT READ CORRECTLY ****" << endl;
       exit(0);  
   }
   else
   {cout << "Data File read correctly" << endl;} 
 END_CALCS   

int header;
!!header=1;

//##################################################################################################
PARAMETER_SECTION //##################################################################################################
//##################################################################################################

 LOCAL_CALCS
  const double Linf_LO=set_Linf(2); const double Linf_HI=set_Linf(3); const double Linf_PH=set_Linf(4);
  const double K_LO=set_K(2); const double K_HI=set_K(3); const double K_PH=set_K(4);
  const double t0_LO=set_t0(2); const double t0_HI=set_t0(3); const double t0_PH=set_t0(4);  
  const double len_cv_LO=set_len_cv(2); const double len_cv_HI=set_len_cv(3); const double len_cv_PH=set_len_cv(4); 
   
  const double Linf_L_LO=set_Linf_L(2); const double Linf_L_HI=set_Linf_L(3); const double Linf_L_PH=set_Linf_L(4);
  const double K_L_LO=set_K_L(2); const double K_L_HI=set_K_L(3); const double K_L_PH=set_K_L(4);
  const double t0_L_LO=set_t0_L(2); const double t0_L_HI=set_t0_L(3); const double t0_L_PH=set_t0_L(4);  
  const double len_cv_L_LO=set_len_cv_L(2); const double len_cv_L_HI=set_len_cv_L(3); const double len_cv_L_PH=set_len_cv_L(4);
  
  const double Linf_F_LO=set_Linf_F(2); const double Linf_F_HI=set_Linf_F(3); const double Linf_F_PH=set_Linf_F(4);
  const double K_F_LO=set_K_F(2); const double K_F_HI=set_K_F(3); const double K_F_PH=set_K_F(4);
  const double t0_F_LO=set_t0_F(2); const double t0_F_HI=set_t0_F(3); const double t0_F_PH=set_t0_F(4);  
  const double len_cv_F_LO=set_len_cv_F(2); const double len_cv_F_HI=set_len_cv_F(3); const double len_cv_F_PH=set_len_cv_F(4);
  
  const double M_constant_LO=set_M_constant(2); const double M_constant_HI=set_M_constant(3); const double M_constant_PH=set_M_constant(4);        
  const double steep_LO=set_steep(2); const double steep_HI=set_steep(3); const double steep_PH=set_steep(4);
  const double log_R0_LO=set_log_R0(2); const double log_R0_HI=set_log_R0(3); const double log_R0_PH=set_log_R0(4);
  const double R_autocorr_LO=set_R_autocorr(2); const double R_autocorr_HI=set_R_autocorr(3); const double R_autocorr_PH=set_R_autocorr(4);
  const double rec_sigma_LO=set_rec_sigma(2); const double rec_sigma_HI=set_rec_sigma(3); const double rec_sigma_PH=set_rec_sigma(4);
  
  const double log_dm_comm_lc_LO=set_log_dm_comm_lc(2); const double log_dm_comm_lc_HI=set_log_dm_comm_lc(3); const double log_dm_comm_lc_PH=set_log_dm_comm_lc(4);
  const double log_dm_GR_ac_LO=set_log_dm_GR_ac(2); const double log_dm_GR_ac_HI=set_log_dm_GR_ac(3); const double log_dm_GR_ac_PH=set_log_dm_GR_ac(4);
  
  const double selpar_A50_comm1_LO=set_selpar_A50_comm1(2); const double selpar_A50_comm1_HI=set_selpar_A50_comm1(3); const double selpar_A50_comm1_PH=set_selpar_A50_comm1(4);
  const double selpar_slope_comm1_LO=set_selpar_slope_comm1(2); const double selpar_slope_comm1_HI=set_selpar_slope_comm1(3); const double selpar_slope_comm1_PH=set_selpar_slope_comm1(4);
  
  const double selpar_A50_GR1_LO=set_selpar_A50_GR1(2); const double selpar_A50_GR1_HI=set_selpar_A50_GR1(3); const double selpar_A50_GR1_PH=set_selpar_A50_GR1(4);
  const double selpar_slope_GR1_LO=set_selpar_slope_GR1(2); const double selpar_slope_GR1_HI=set_selpar_slope_GR1(3); const double selpar_slope_GR1_PH=set_selpar_slope_GR1(4);
  const double selpar_A50_GR2_LO=set_selpar_A50_GR2(2); const double selpar_A50_GR2_HI=set_selpar_A50_GR2(3); const double selpar_A50_GR2_PH=set_selpar_A50_GR2(4);
  const double selpar_slope_GR2_LO=set_selpar_slope_GR2(2); const double selpar_slope_GR2_HI=set_selpar_slope_GR2(3); const double selpar_slope_GR2_PH=set_selpar_slope_GR2(4);  
  
  const double log_q_HB_LO=set_log_q_HB(2); const double log_q_HB_HI=set_log_q_HB(3); const double log_q_HB_PH=set_log_q_HB(4);
  
  const double F_init_LO=set_F_init(2); const double F_init_HI=set_F_init(3); const double F_init_PH=set_F_init(4);
  const double log_avg_F_comm_LO=set_log_avg_F_comm(2); const double log_avg_F_comm_HI=set_log_avg_F_comm(3); const double log_avg_F_comm_PH=set_log_avg_F_comm(4);
  const double log_avg_F_GR_LO=set_log_avg_F_GR(2); const double log_avg_F_GR_HI=set_log_avg_F_GR(3); const double log_avg_F_GR_PH=set_log_avg_F_GR(4); 
  
  //-dev vectors-----------------------------------------------------------------------------------------------------------  
  const double log_F_dev_comm_LO=set_log_F_dev_comm(1); const double log_F_dev_comm_HI=set_log_F_dev_comm(2); const double log_F_dev_comm_PH=set_log_F_dev_comm(3);  
  const double log_F_dev_GR_LO=set_log_F_dev_GR(1); const double log_F_dev_GR_HI=set_log_F_dev_GR(2); const double log_F_dev_GR_PH=set_log_F_dev_GR(3);     
  
  const double log_RWq_LO=set_log_RWq_dev(1); const double log_RWq_HI=set_log_RWq_dev(2); const double log_RWq_PH=set_log_RWq_dev(3);  
  
  const double log_rec_dev_LO=set_log_rec_dev(1); const double log_rec_dev_HI=set_log_rec_dev(2); const double log_rec_dev_PH=set_log_rec_dev(3);          
  const double log_Nage_dev_LO=set_log_Nage_dev(1); const double log_Nage_dev_HI=set_log_Nage_dev(2); const double log_Nage_dev_PH=set_log_Nage_dev(3);          
  
 END_CALCS
 
////--------------Growth--------------------------------------------------------------------------- 
  //Population growth parms and conversions
  init_bounded_number Linf(Linf_LO,Linf_HI,Linf_PH);
  init_bounded_number K(K_LO,K_HI,K_PH);
  init_bounded_number t0(t0_LO,t0_HI,t0_PH);
  init_bounded_number len_cv_val(len_cv_LO,len_cv_HI,len_cv_PH);  
  vector Linf_out(1,8);
  vector K_out(1,8);
  vector t0_out(1,8);
  vector len_cv_val_out(1,8);

  vector meanlen_TL(1,nages);   //mean total length (mm) at age all fish
   
  vector wgt_g(1,nages);        //whole wgt in g
  vector wgt_kg(1,nages);       //whole wgt in kg
  vector wgt_mt(1,nages);       //whole wgt in mt
  vector wgt_klb(1,nages);      //whole wgt in 1000 lb
  vector wgt_lb(1,nages);       //whole wgt in lb  
  
  init_bounded_number Linf_L(Linf_L_LO,Linf_L_HI,Linf_L_PH);
  init_bounded_number K_L(K_L_LO,K_L_HI,K_L_PH);
  init_bounded_number t0_L(t0_L_LO,t0_L_HI,t0_L_PH);
  init_bounded_number len_cv_val_L(len_cv_L_LO,len_cv_L_HI,len_cv_L_PH);  
  vector Linf_L_out(1,8);
  vector K_L_out(1,8);
  vector t0_L_out(1,8);
  vector len_cv_val_L_out(1,8);
  vector meanlen_TL_L(1,nages);   //mean total length (mm) at age all fish
  
  vector wgt_g_L(1,nages);        //whole wgt in g
  vector wgt_kg_L(1,nages);       //whole wgt in kg
  vector wgt_mt_L(1,nages);       //whole wgt in mt
  vector wgt_klb_L(1,nages);      //whole wgt in 1000 lb
  vector wgt_lb_L(1,nages);       //whole wgt in lb  
  vector wgt_klb_gut_L(1,nages);  //gutted wgt in 1000 lb  
  vector wgt_lb_gut_L(1,nages);   //gutted wgt in lb
     
  init_bounded_number Linf_F(Linf_F_LO,Linf_F_HI,Linf_F_PH);
  init_bounded_number K_F(K_F_LO,K_F_HI,K_F_PH);
  init_bounded_number t0_F(t0_F_LO,t0_F_HI,t0_F_PH);
  init_bounded_number len_cv_val_F(len_cv_F_LO,len_cv_F_HI,len_cv_F_PH);  
  vector Linf_F_out(1,8);
  vector K_F_out(1,8);
  vector t0_F_out(1,8);
  vector len_cv_val_F_out(1,8);
  vector meanlen_TL_F(1,nages);   //mean total length (mm) at age all fish
  
  vector wgt_g_F(1,nages);        //whole wgt in g
  vector wgt_kg_F(1,nages);       //whole wgt in kg
  vector wgt_mt_F(1,nages);       //whole wgt in mt
  vector wgt_klb_F(1,nages);      //whole wgt in 1000 lb
  vector wgt_lb_F(1,nages);       //whole wgt in lb  

  
  //vector batchfec(1,nages);     //batch fecundity at age	
  //vector fec(1,nages);          //annual fecundity at age	
  
  matrix len_comm_mm(styr,endyr,1,nages);          //mean length at age of commercial handline landings in mm 
  matrix wholewgt_comm_klb(styr,endyr,1,nages);    //whole wgt of commercial handline landings in 1000 lb   
  //matrix len_cL_mm(styr,endyr,1,nages);          //mean length at age of commercial longline landings in mm 
  //matrix wholewgt_cL_klb(styr,endyr,1,nages);    //whole wgt of commercial longline landings in 1000 lb   
  matrix len_HB_mm(styr,endyr,1,nages);          //mean length at age of HB landings in mm     
  matrix wholewgt_HB_klb(styr,endyr,1,nages);    //whole wgt of HB landings in 1000 lb  
  matrix len_GR_mm(styr,endyr,1,nages);          //mean length at age of GR landings in mm     
  matrix wholewgt_GR_klb(styr,endyr,1,nages);    //whole wgt of GR landings in 1000 lb    
  
  matrix lenprob(1,nages,1,nlenbins);           //distn of size at age (age-length key, 3 cm bins) in population
  number zscore_len;                            //standardized normal values used for computing lenprob
  vector cprob_lenvec(1,nlenbins);              //cumulative probabilities used for computing lenprob
  number zscore_lzero;                          //standardized normal values for length = 0
  number cprob_lzero;                           //length probability mass below zero, used for computing lenprob
   
  matrix lenprob_L(1,nages,1,nlenbins); 
  number zscore_len_L;                            //standardized normal values used for computing lenprob
  vector cprob_lenvec_L(1,nlenbins);              //cumulative probabilities used for computing lenprob
  number zscore_lzero_L;                          //standardized normal values for length = 0
  number cprob_lzero_L;                           //length probability mass below zero, used for computing lenprob
  
  matrix lenprob_F(1,nages,1,nlenbins);
  number zscore_len_F;                            //standardized normal values used for computing lenprob
  vector cprob_lenvec_F(1,nlenbins);              //cumulative probabilities used for computing lenprob
  number zscore_lzero_F;                          //standardized normal values for length = 0
  number cprob_lzero_F;                           //length probability mass below zero, used for computing lenprob
  
  
  //matrices below are used to match length comps
  matrix lenprob_comm(1,nages,1,nlenbins);     //distn of size at age in comm
  //matrix lenprob_cL(1,nages,1,nlenbins);     //distn of size at age in cL
  matrix lenprob_HB(1,nages,1,nlenbins);     //distn of size at age in HB  
  matrix lenprob_GR(1,nages,1,nlenbins);     //distn of size at age in GR  
  
  vector len_sd(1,nages);
  vector len_cv(1,nages); //for fishgraph 
  //All Fishery-dependent
  vector len_sd_L(1,nages);
  vector len_cv_L(1,nages); //for fishgraph 
  //Females
  vector len_sd_F(1,nages);
  vector len_cv_F(1,nages);
  
//----Predicted length and age compositions
  matrix pred_comm_lenc(1,nyr_comm_lenc,1,nlenbins); //predicted length comps pooled across years
  matrix pred_comm_lenc_yr(1,nyr_comm_lenc_pool,1,nlenbins); //annual predicted length comps
  //matrix pred_cL_lenc(1,nyr_cL_lenc,1,nlenbins);
  //matrix pred_HB_lenc(1,nyr_HB_lenc,1,nlenbins);   
  //matrix pred_GR_lenc(1,nyr_GR_lenc,1,nlenbins);   
  matrix pred_GR_agec(1,nyr_GR_agec,1,nages_agec);  
  matrix pred_GR_agec_allages(1,nyr_GR_agec,1,nages);  
  matrix ErrorFree_GR_agec(1,nyr_GR_agec,1,nages);  
  
//Sample size (perhaps adjusted herein) used in fitting comp data
  vector nsamp_comm_lenc_allyr(styr,endyr);
  //vector nsamp_cL_lenc_allyr(styr,endyr);
  // vector nsamp_HB_lenc_allyr(styr,endyr);  
  //vector nsamp_GR_lenc_allyr(styr,endyr);
  vector nsamp_GR_agec_allyr(styr,endyr);  
  
//Nfish used in MCB analysis (not used in fitting)
  vector nfish_comm_lenc_allyr(styr,endyr);
  //vector nfish_cL_lenc_allyr(styr,endyr);
  // vector nfish_HB_lenc_allyr(styr,endyr);
  //vector nfish_GR_lenc_allyr(styr,endyr);
  vector nfish_GR_agec_allyr(styr,endyr);  

//Computed effective sample size for output (not used in fitting)
  vector neff_comm_lenc_allyr(styr,endyr);
  //vector neff_cL_lenc_allyr(styr,endyr);
  // vector neff_HB_lenc_allyr(styr,endyr); 
  //vector neff_GR_lenc_allyr(styr,endyr);
  vector neff_GR_agec_allyr(styr,endyr);  

//-----Population-----------------------------------------------------------------------------------
  matrix N(styr,endyr+1,1,nages);           //Population numbers by year and age at start of yr
  matrix N_mdyr(styr,endyr,1,nages);        //Population numbers by year and age at mdpt of yr: used for comps and cpue
  matrix N_spawn(styr,endyr,1,nages);       //Population numbers by year and age at peaking spawning: used for SSB  
  init_bounded_vector log_Nage_dev(2,nages,log_Nage_dev_LO,log_Nage_dev_HI,log_Nage_dev_PH);
  vector log_Nage_dev_output(1,nages);      //used in output. equals zero for first age
  matrix B(styr,endyr+1,1,nages);           //Population biomass by year and age at start of yr
  vector totB(styr,endyr+1);                //Total biomass by year
  vector totN(styr,endyr+1);                //Total abundance by year
  vector SSB(styr,endyr);                   //Total spawning biomass by year (female mature biomass)
  vector SSB_knum(styr,endyr);              //Total spawning numbers by year (number of mature females)  
  vector rec(styr,endyr+1);                 //Recruits by year
  vector prop_f(1,nages);
  //vector prop_m(1,nages);
  vector maturity_f(1,nages);
  //vector maturity_m(1,nages);
  vector reprod(1,nages);
  vector reprodknum(1,nages);
 
//---Stock-Recruit Function (Beverton-Holt, steepness parameterization)----------
  init_bounded_number log_R0(log_R0_LO,log_R0_HI,log_R0_PH);        //log(virgin Recruitment)
  vector log_R0_out(1,8);
  number R0;                                  //virgin recruitment
  init_bounded_number steep(steep_LO,steep_HI,steep_PH); //steepness
  vector steep_out(1,8);
  init_bounded_number rec_sigma(rec_sigma_LO,rec_sigma_HI,rec_sigma_PH);  //sd recruitment residuals  
  vector rec_sigma_out(1,8);
  init_bounded_number R_autocorr(R_autocorr_LO,R_autocorr_HI,R_autocorr_PH);  //autocorrelation in SR  
  vector R_autocorr_out(1,8);

  number rec_sigma_sq;                        //square of rec_sigma      
  number rec_logL_add;                        //additive term in -logL term   
 
  init_bounded_dev_vector log_rec_dev(styr_rec_dev,endyr_rec_dev,log_rec_dev_LO,log_rec_dev_HI,log_rec_dev_PH);
  vector log_rec_dev_output(styr,endyr+1);             //used in t.series output. equals zero except for yrs in log_rec_dev
  vector log_rec_dev_out(styr_rec_dev,endyr_rec_dev);  //used in output for bound checking
  
  number var_rec_dev;                           //variance of log recruitment deviations, from yrs with unconstrainted S-R(XXXX-XXXX)
  number sigma_rec_dev;                         //sample SD of log residuals (may not equal rec_sigma 
  number BiasCor;                               //Bias correction in equilibrium recruits
  number S0;                                    //equal to spr_F0*R0 = virgin SSB
  number B0;                                    //equal to bpr_F0*R0 = virgin B  
  sdreport_number R1;                                    //Recruits in styr
  number R_virgin;                              //unfished recruitment with bias correction
  vector SdS0(styr,endyr);                      //Spawners relative to the unfished level              

  
  init_bounded_number log_dm_comm_lc(log_dm_comm_lc_LO,log_dm_comm_lc_HI,log_dm_comm_lc_PH);
  init_bounded_number log_dm_GR_ac(log_dm_GR_ac_LO,log_dm_GR_ac_HI,log_dm_GR_ac_PH);

  vector log_dm_comm_lc_out(1,8);
  vector log_dm_GR_ac_out(1,8);
//-----------------------------------------------------------------------------------------------------------------------------------------------
////---Selectivity-------------------------------------------------------------------------

//Commercial handline-------------------------------------------------
  matrix sel_comm(styr,endyr,1,nages);
  vector sel_comm_vec(1,nages);
  
  init_bounded_number selpar_A50_comm1(selpar_A50_comm1_LO,selpar_A50_comm1_HI,selpar_A50_comm1_PH);
  init_bounded_number selpar_slope_comm1(selpar_slope_comm1_LO,selpar_slope_comm1_HI,selpar_slope_comm1_PH);
  
  vector selpar_A50_comm1_out(1,8);
  vector selpar_slope_comm1_out(1,8);

  //Headboat -------------------------------------------------
  matrix sel_HB(styr,endyr,1,nages); // Still need to define sel_HB to associate with HB index, but can just set equal to sel_GR below
  vector sel_HB_block1(1,nages);
  vector sel_HB_block2(1,nages);
  
  //General Rec 
  matrix sel_GR(styr,endyr,1,nages);
  vector sel_GR_block1(1,nages);
  vector sel_GR_block2(1,nages);
 
  init_bounded_number selpar_A50_GR1(selpar_A50_GR1_LO,selpar_A50_GR1_HI,selpar_A50_GR1_PH);
  init_bounded_number selpar_slope_GR1(selpar_slope_GR1_LO,selpar_slope_GR1_HI,selpar_slope_GR1_PH);
  init_bounded_number selpar_A50_GR2(selpar_A50_GR2_LO,selpar_A50_GR2_HI,selpar_A50_GR2_PH);
  init_bounded_number selpar_slope_GR2(selpar_slope_GR2_LO,selpar_slope_GR2_HI,selpar_slope_GR2_PH);
  
  vector selpar_A50_GR1_out(1,8);
  vector selpar_slope_GR1_out(1,8);
  vector selpar_A50_GR2_out(1,8);
  vector selpar_slope_GR2_out(1,8);

//Weighted total selectivity--------------------------------------------  
  //effort-weighted, recent selectivities
  vector sel_wgted_L(1,nages);  //toward landings 
  vector sel_wgted_tot(1,nages);//toward Z, landings plus deads discards

//-----------------------------------------------------------------------------------------------------------------------------------------------
//-------CPUE Predictions--------------------------------
  //vector pred_comm_cpue(styr_comm_cpue,endyr_comm_cpue);                   //predicted comm index (weight fish per effort)
  //matrix N_comm(styr_comm_cpue,endyr_comm_cpue,1,nages);                   //used to compute comm index
  //vector pred_cL_cpue(styr_cL_cpue,endyr_cL_cpue);                   //predicted cL index (weight fish per effort)
 // matrix N_cL(styr_cL_cpue,endyr_cL_cpue,1,nages);                   //used to compute cL index
  vector pred_HB_cpue(styr_HB_cpue,endyr_HB_cpue);                   //predicted HB index (number fish per effort)
  matrix N_HB(styr_HB_cpue,endyr_HB_cpue,1,nages);                   //used to compute HB index
  
//---Catchability (CPUE q's)----------------------------------------------------------
  //init_bounded_number log_q_comm(log_q_comm_LO,log_q_comm_HI,log_q_comm_PH);
  //init_bounded_number log_q_cL(log_q_cL_LO,log_q_cL_HI,log_q_cL_PH);  
  init_bounded_number log_q_HB(log_q_HB_LO,log_q_HB_HI,log_q_HB_PH);
  
  //vector log_q_comm_out(1,8);
 // vector log_q_cL_out(1,8);  
  vector log_q_HB_out(1,8);
  
  number q_rate;
  //vector q_rate_fcn_comm(styr_comm_cpue,endyr_comm_cpue);         //increase due to technology creep (saturates in 2003)
  //vector q_rate_fcn_cL(styr_cL_cpue,endyr_cL_cpue);         //increase due to technology creep (saturates in 2003)  
  vector q_rate_fcn_HB(styr_HB_cpue,endyr_HB_cpue);         //increase due to technology creep (saturates in 2003) 
  
//  init_bounded_number q_DD_beta(0.1,0.9,set_q_DD_phase);    //not estimated so commented out and declared as number (below)
  number q_DD_beta;
  vector q_DD_fcn(styr,endyr);    //density dependent function as a multiple of q (scaled a la Katsukawa and Matsuda. 2003)
  number B0_q_DD;                 //B0 of ages q_DD_age plus
  vector B_q_DD(styr,endyr);      //annual biomass of ages q_DD_age plus

//Fishery dependent random walk catchability
 //init_bounded_vector q_RW_log_dev_comm(styr_comm_cpue,endyr_comm_cpue-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
 //init_bounded_vector q_RW_log_dev_cL(styr_cL_cpue,endyr_cL_cpue-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
 init_bounded_vector q_RW_log_dev_HB(styr_HB_cpue,endyr_HB_cpue-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);

//Fishery dependent catchability over time, may be constant
 //vector q_comm(styr_comm_cpue,endyr_comm_cpue); 
 //vector q_cL(styr_cL_cpue,endyr_cL_cpue);
 vector q_HB(styr_HB_cpue,endyr_HB_cpue); 

//----------------------------------------------------------------------------------------------------------------------------------------------- 
//---Landings in numbers (total or 1000 fish) and in wgt (whole klb)--------------------------------------------------
  matrix L_comm_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_comm_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age    
  vector pred_comm_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages  
  vector pred_comm_L_klb(styr,endyr);      //yearly landings in 1000 lb whole summed over ages

  //matrix L_cL_num(styr,endyr,1,nages);   //landings (numbers) at age
  //matrix L_cL_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age    
  //vector pred_cL_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages  
  //vector pred_cL_L_klb(styr,endyr);      //yearly landings in 1000 lb whole summed over ages

  matrix L_GR_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_GR_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age
  vector pred_GR_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages
  vector pred_GR_L_klb(styr,endyr);      //yearly landings in 1000 lb whole summed over ages

  matrix L_total_num(styr,endyr,1,nages);//total landings in number at age
  matrix L_total_klb(styr,endyr,1,nages);//landings in klb whole wgt at age 
  vector L_total_knum_yr(styr,endyr);    //total landings in 1000 fish by yr summed over ages  
  vector L_total_klb_yr(styr,endyr);     //total landings (klb whole wgt) by yr summed over ages

////---MSY calcs----------------------------------------------------------------------------
  number F_comm_prop;       //proportion of F_sum attributable to comm, last X=selpar_n_yrs_wgted yrs
  //number F_cL_prop;       //proportion of F_sum attributable to comm, last X=selpar_n_yrs_wgted yrs  
  number F_GR_prop;       //proportion of F_sum attributable to GR, last X=selpar_n_yrs_wgted yrs
  
  number F_init_comm_prop;  //proportion of F_init attributable to comm, first X yrs
  //number F_init_cL_prop;  //proportion of F_init attributable to cL, first X yrs
  number F_init_GR_prop;  //proportion of F_init attributable to GR, first X yrs
  
  number F_temp_sum;      //sum of geom mean Fsum's in last X yrs, used to compute F_fishery_prop

  vector F_end(1,nages);
  vector F_end_L(1,nages);     
  number F_end_apex;
  
  number SSB_msy_out;           //SSB (total mature biomass) at msy
  number F_msy_out;             //F at msy
  number msy_klb_out;           //max sustainable yield (1000 lb whole wgt)
  number msy_knum_out;          //max sustainable yield (1000 fish)  
  number B_msy_out;             //total biomass at MSY 
  number R_msy_out;             //equilibrium recruitment at F=Fmsy
  number spr_msy_out;           //spr at F=Fmsy

  number F20_dum;				//intermediate calculation for F20
  number F30_dum;				//intermediate calculation for F30
  number F40_dum;				//intermediate calculation for F40
  number F20_out;              	//F20
  number F30_out;              	//F30
  number F40_out;              	//F40
  number SSB_F30_out;
  number SSB_F30_knum_out;  
  number B_F30_out;
  number R_F30_out;
  number L_F30_knum_out;
  number L_F30_klb_out; 
 
  number SSB_F40_out; 
  number SSB_F40_knum_out;    
  number B_F40_out;
  number R_F40_out;
  number L_F40_knum_out;
  number L_F40_klb_out;   
  number rec_mean;  			//arithmetic average recruitment used in SPR-related quantities

  vector N_age_msy(1,nages);         //numbers at age for MSY calculations: beginning of yr
  vector N_age_msy_spawn(1,nages);   //numbers at age for MSY calculations: time of peak spawning  
  vector L_age_msy(1,nages);         //landings at age for MSY calculations
  vector Z_age_msy(1,nages);         //total mortality at age for MSY calculations
  vector F_L_age_msy(1,nages);       //fishing mortality landings (not discards) at age for MSY calculations
  vector F_msy(1,n_iter_msy);        //values of full F to be used in equilibrium calculations
  vector spr_msy(1,n_iter_msy);      //reproductive capacity-per-recruit values corresponding to F values in F_msy
  vector R_eq(1,n_iter_msy);         //equilibrium recruitment values corresponding to F values in F_msy
  vector L_eq_klb(1,n_iter_msy);     //equilibrium landings(klb whole wgt) values corresponding to F values in F_msy
  vector L_eq_knum(1,n_iter_msy);    //equilibrium landings(1000 fish) values corresponding to F values in F_msy
  vector SSB_eq(1,n_iter_msy);       //equilibrium reproductive capacity values corresponding to F values in F_msy
  vector SSB_eq_knum(1,n_iter_msy);
  vector B_eq(1,n_iter_msy);         //equilibrium biomass values corresponding to F values in F_msy
  
  vector FdF_msy(styr,endyr);
  vector FdF30(styr,endyr);
  vector FdF40(styr,endyr);
  vector SdSSB_msy(styr,endyr);	 
  number SdSSB_msy_end;
  number FdF_msy_end;
  number FdF_msy_end_mean;           //geometric mean of last X yrs  
  
  vector SdSSB_F30(styr,endyr);	 
  vector Sdmsst_F30(styr,endyr);	 
  number SdSSB_F30_end;
  number Sdmsst_F30_end;
  number FdF30_end_mean;             //geometric mean of last selpar_n_yrs_wgted yrs  
  vector L_age_F30(1,nages);         //landings at age for F30 calculations
  
  vector SdSSB_F40(styr,endyr);	 
  vector Sdmsst_F40(styr,endyr);	 
  number SdSSB_F40_end;
  number Sdmsst_F40_end;
  number FdF40_end_mean;             //geometric mean of last selpar_n_yrs_wgted yrs  
  number Fend_mean_temp;			 //intermediate calc for geometric mean of last selpar_n_yrs_wgted yrs
  number Fend_mean;					 //geometric mean of last selpar_n_yrs_wgted yrs
  vector L_age_F40(1,nages);         //landings at age for F40 calculations
  
  vector wgt_wgted_L_klb(1,nages);   //fishery-weighted average weight at age of landings in whole weight 
  number wgt_wgted_L_denom;          //used in intermediate calculations

  number iter_inc_msy;               //increments used to compute msy, equals 1/(n_iter_msy-1)
  
////--------Mortality------------------------------------------------------------------

  vector M(1,nages);                         //age-dependent natural mortality
  init_bounded_number M_constant(M_constant_LO,M_constant_HI,M_constant_PH);   //age-indpendent: used only for MSST
  vector M_constant_out(1,8);
  number smsy2msstM;                         //scales Smsy to get msst using (1-M). Used only in output.
  number smsy2msst75;                        //scales Smsy to get msst using 75%. Used only in output.  
  
  matrix F(styr,endyr,1,nages);
  vector Fsum(styr,endyr);                   //Full fishing mortality rate by year
  vector Fapex(styr,endyr);                  //Max across ages, fishing mortality rate by year (may differ from Fsum bc of dome-shaped sel 
  matrix Z(styr,endyr,1,nages);

  init_bounded_number log_avg_F_comm(log_avg_F_comm_LO,log_avg_F_comm_HI,log_avg_F_comm_PH);
  vector log_avg_F_comm_out(1,8);  
  init_bounded_dev_vector log_F_dev_comm(styr_comm_L,endyr_comm_L,log_F_dev_comm_LO,log_F_dev_comm_HI,log_F_dev_comm_PH);
  vector log_F_dev_comm_out(styr_comm_L,endyr_comm_L);
  matrix F_comm(styr,endyr,1,nages);
  vector F_comm_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_comm;
  number log_F_dev_end_comm; 

  init_bounded_number log_avg_F_GR(log_avg_F_GR_LO,log_avg_F_GR_HI,log_avg_F_GR_PH);
  vector log_avg_F_GR_out(1,8); 
  init_bounded_dev_vector log_F_dev_GR(styr_GR_L,endyr_GR_L,log_F_dev_GR_LO,log_F_dev_GR_HI,log_F_dev_GR_PH);    
  vector log_F_dev_GR_out(styr_GR_L,endyr_GR_L);
  matrix F_GR(styr,endyr,1,nages);
  vector F_GR_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_GR;    
  number log_F_dev_end_GR;

  init_bounded_number F_init(F_init_LO,F_init_HI,F_init_PH); //scales early F for initialization
  vector F_init_out(1,8); 
  number F_init_denom;  //interim calculation. From Erik's red snapper ASPM

  //number F_init_ratio;   //scales initial F, which is read in as a fixed value
  //vector sel_initial(1,nages);       //initial selectivity (a combination of recreational and commercial selectivities)  

//---Per-recruit stuff----------------------------------------------------------------------------------
  vector N_age_spr(1,nages);         //numbers at age for SPR calculations: beginning of year
  vector N_age_spr_spawn(1,nages);   //numbers at age for SPR calculations: time of peak spawning  
  vector L_age_spr(1,nages);         //catch at age for SPR calculations
  vector Z_age_spr(1,nages);         //total mortality at age for SPR calculations
  vector spr_static(styr,endyr);     //vector of static SPR values by year
  vector F_L_age_spr(1,nages);       //fishing mortality of landings (not discards) at age for SPR calculations
  vector F_spr(1,n_iter_spr);        //values of full F to be used in per-recruit calculations
  vector spr_spr(1,n_iter_spr);      //reproductive capacity-per-recruit values corresponding to F values in F_spr
  vector spr_ratio(1,n_iter_spr);    //reproductive capacity-per-recruit relative to spr_F0 values corresponding to F values in F_spr
  vector L_spr(1,n_iter_spr);        //landings(lb)-per-recruit (ypr) values corresponding to F values in F_spr

  vector N_spr_F0(1,nages);          //Used to compute spr at F=0: at time of peak spawning
  vector N_bpr_F0(1,nages);          //Used to compute bpr at F=0: at start of year  
  vector N_spr_initial(1,nages);     //Initial spawners per recruit at age given initial F
  vector N_initial_eq(1,nages);      //Initial equilibrium abundance at age
  vector F_initial(1,nages);         //initial F at age
  vector Z_initial(1,nages);         //initial Z at age
  number spr_initial;                //initial spawners per recruit
  number spr_F0;                     //Spawning biomass per recruit at F=0
  number bpr_F0;                     //Biomass per recruit at F=0

  number iter_inc_spr;               //increments used to compute msy, equals max_F_spr_msy/(n_iter_spr-1)

////-------SDNR output-----------------------------------------------------------------------------
 
  number sdnr_lc_comm;
  //number sdnr_lc_cL;
  //number sdnr_lc_HB; 
  //number sdnr_lc_GR;   
  number sdnr_ac_GR;
 // number sdnr_I_comm;
 // number sdnr_I_cL;  
  number sdnr_I_HB; 
    
////-------Objective function components-----------------------------------------------------------------------------
  number w_L;
  
 // number w_I_comm;
 // number w_I_cL;  
  number w_I_HB;
   
  number w_lc_comm;
  //number w_lc_cL; 
  //number w_lc_HB;
  //number w_lc_GR; 
  number w_ac_GR;  
  
  number w_Nage_init;  
  number w_rec;
  number w_rec_early;
  number w_rec_end;
  number w_fullF;  
  number w_Ftune;

  number f_comm_L; 
  //number f_cL_L; 
  // number f_HB_L; 
  number f_GR_L; 

  //number f_comm_cpue;
  //number f_cL_cpue;  
  number f_HB_cpue;
 
  //number f_comm_RWq_cpue;
  //number f_cL_RWq_cpue;
  number f_HB_RWq_cpue;
  
  number f_comm_lenc;
  //number f_cL_lenc;  
  // number f_HB_lenc;
  number f_GR_lenc;

  number f_GR_agec;

//Likelihood Components for Data (for WAIC Calcs)
  vector LL_HB_cpue(styr_HB_cpue,endyr_HB_cpue)
  vector LL_comm_L(styr_comm_L,endyr_comm_L)
  vector LL_GR_L(styr_GR_L,endyr_GR_L)
  number LL_comm_lc
  vector LL_GR_ac(1,26)
  
//  Penalties and constraints. Not all are used.
  number f_Nage_init;              //weight on log devs to estimate initial abundance (excluding first age)
  number f_rec_dev;                //weight on recruitment deviations to fit S-R curve
  number f_rec_dev_early;          //extra weight on deviations in first recruitment stanza
  number f_rec_dev_end;            //extra weight on deviations in ending recruitment stanza
  number f_fullF_constraint;       //penalty for Fapex>X
  number f_Ftune;                  //penalty for tuning F in Ftune yr.  Not applied in final optimization phase.
  number f_priors;                 //prior information on parameters
  
   //init_number xdum;
   objective_function_value fval;
   number fval_data;
   number grad_max;
  
 //--Dummy variables ----
   number denom;                   //denominator used in some calculations
   number numer;                   //numerator used in some calculations

  //---------- Projection quantities--------------------------------------------------------------
  number F_reg_proj;						   //value used to define the projections				
  vector F_proj(styr_proj,endyr_proj);         //F by yr for projections (=F_reg_proj after regulations start, current F till then)  
  vector L_knum_proj(styr_proj,endyr_proj);    //total landings in 1000 fish for projections  
  vector L_klb_proj(styr_proj,endyr_proj);     //total landings in weight (1000 lb) for projections

  vector B_proj(styr_proj,endyr_proj);         //Biomass for projections  
  vector SSB_proj(styr_proj,endyr_proj);       //SSB for projections  
  vector R_proj(styr_proj,endyr_proj);     	   //recruits for projections
  vector FL_age_proj(1,nages);      		   //F (landings) by age for projections       
  
  matrix N_proj(styr_proj,endyr_proj,1,nages);           //Population numbers by year and age at start of yr
  matrix N_spawn_proj(styr_proj,endyr_proj,1,nages);     //Population numbers by year and age at peaking spawning: used for SSB in projections 
  matrix Z_proj(styr_proj,endyr_proj,1,nages);           //Z by year and age for projections 
  matrix L_age_proj(styr_proj,endyr_proj,1,nages);       //Projected landings at age in numbers    
  
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//INITIALIZATION_SECTION
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
GLOBALS_SECTION
  #include "admodel.h"          // Include AD class definitions
  #include "admb2r.cpp"         // Include S-compatible output functions (needs preceding)
  #include <time.h>
	time_t start,finish;
	long hour,minute,second;	
	double elapsed_time;
	
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
RUNTIME_SECTION
 maximum_function_evaluations 1000, 2000,3000, 5000, 10000;//, 10000, 10000;
 convergence_criteria 1e-2, 1e-2,1e-3, 1e-3, 1e-4;//, 1e-4, 1e-4;
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PRELIMINARY_CALCS_SECTION

// Set values of fixed parameters or set initial guess of estimated parameters
 
 //Population		
  Linf=set_Linf(1);
  K=set_K(1);
  t0=set_t0(1);
  len_cv_val=set_len_cv(1);
  
  //All fisheries
  Linf_L=set_Linf_L(1);  
  K_L=set_K_L(1);
  t0_L=set_t0_L(1);
  len_cv_val_L=set_len_cv_L(1);
  //Females
  Linf_F=set_Linf_F(1);  
  K_F=set_K_F(1);
  t0_F=set_t0_F(1);
  len_cv_val_F=set_len_cv_F(1);
  
  M=set_M; 
  M_constant=set_M_constant(1);
  smsy2msstM=1.0-M_constant;
  smsy2msst75=0.75;  
  
  log_R0=set_log_R0(1);
  steep=set_steep(1);
  R_autocorr=set_R_autocorr(1);
  rec_sigma=set_rec_sigma(1);
  
  log_dm_comm_lc=set_log_dm_comm_lc(1);
  //log_dm_cL_lc=set_log_dm_cL_lc(1);
  // log_dm_HB_lc=set_log_dm_HB_lc(1);
  //log_dm_GR_lc=set_log_dm_GR_lc(1);
  log_dm_GR_ac=set_log_dm_GR_ac(1);
  
 // log_q_comm=set_log_q_comm(1);
  //log_q_cL=set_log_q_cL(1);  
  log_q_HB=set_log_q_HB(1);
  
  q_rate=set_q_rate;
  //q_rate_fcn_comm=1.0;
  //q_rate_fcn_cL=1.0;    
  q_rate_fcn_HB=1.0;   
  q_DD_beta=set_q_DD_beta;
  q_DD_fcn=1.0;

  //q_RW_log_dev_comm.initialize();
 // q_RW_log_dev_cL.initialize();   
  q_RW_log_dev_HB.initialize(); 
  
   if (set_q_rate_phase<0 & q_rate!=0.0)
  {
      
    for (iyear=styr_HB_cpue; iyear<=endyr_HB_cpue; iyear++)
      {   if (iyear>styr_HB_cpue & iyear <=2003) 
          {//q_rate_fcn_HB(iyear)=(1.0+q_rate)*q_rate_fcn_HB(iyear-1); //compound
             q_rate_fcn_HB(iyear)=(1.0+(iyear-styr_HB_cpue)*q_rate)*q_rate_fcn_HB(styr_HB_cpue);  //linear
          }
          if (iyear>2003) {q_rate_fcn_HB(iyear)=q_rate_fcn_HB(iyear-1);} 
      }   
  } //end q_rate conditional      

  w_L=set_w_L;
  
 // w_I_comm=set_w_I_comm;
 // w_I_cL=set_w_I_cL;
  w_I_HB=set_w_I_HB;
  
  w_lc_comm=set_w_lc_comm;
  //w_lc_cL=set_w_lc_cL;  
  //w_lc_HB=set_w_lc_HB;
  //w_lc_GR=set_w_lc_GR;  
  w_ac_GR=set_w_ac_GR; 
  
  w_Nage_init=set_w_Nage_init;
  w_rec=set_w_rec;
  w_rec_early=set_w_rec_early;
  w_rec_end=set_w_rec_end;
  w_fullF=set_w_fullF;
  w_Ftune=set_w_Ftune;
  
  F_init=set_F_init(1);  

  log_avg_F_comm=set_log_avg_F_comm(1);
  //log_avg_F_cL=set_log_avg_F_cL(1);  
  // log_avg_F_HB=set_log_avg_F_HB(1); 
  log_avg_F_GR=set_log_avg_F_GR(1); 
    
  log_F_dev_comm=set_log_F_dev_comm_vals;
  //log_F_dev_cL=set_log_F_dev_cL_vals;
  // log_F_dev_HB=set_log_F_dev_HB_vals;
  log_F_dev_GR=set_log_F_dev_GR_vals;
 
  selpar_A50_comm1=set_selpar_A50_comm1(1);
  selpar_slope_comm1=set_selpar_slope_comm1(1);
  //selpar_A50_comm2=set_selpar_A50_comm2(1);
  //selpar_slope_comm2=set_selpar_slope_comm2(1);
  // selpar_A502_comm2=set_selpar_A502_comm2(1);
  // selpar_slope2_comm2=set_selpar_slope2_comm2(1);  
  //selpar_A50_comm3=set_selpar_A50_comm3(1);
  //selpar_slope_comm3=set_selpar_slope_comm3(1);  
  
  selpar_A50_GR1=set_selpar_A50_GR1(1);
  selpar_slope_GR1=set_selpar_slope_GR1(1);
  selpar_A50_GR2=set_selpar_A50_GR2(1);
  selpar_slope_GR2=set_selpar_slope_GR2(1);
  // selpar_A502_GR2=set_selpar_A502_GR2(1);
  // selpar_slope2_GR2=set_selpar_slope2_GR2(1);  
  //selpar_A50_GR3=set_selpar_A50_GR3(1);
  //selpar_slope_GR3=set_selpar_slope_GR3(1);  

 sqrt2pi=sqrt(2.*3.14159265);
 g2mt=0.000001;         //conversion of grams to metric tons
 g2kg=0.001;            //conversion of grams to kg 
 mt2klb=2.20462;        //conversion of metric tons to 1000 lb 
 mt2lb=mt2klb*1000.0;   //conversion of metric tons to lb
 g2klb=g2mt*mt2klb;     //conversion of grams to 1000 lb 
 dzero=0.00001;         
 huge_number=1.0e+10;   
 
 SSB_msy_out=0.0;

 iter_inc_msy=max_F_spr_msy/(n_iter_msy-1);
 iter_inc_spr=max_F_spr_msy/(n_iter_spr-1); 

 maturity_f=maturity_f_obs;
 
 prop_f=prop_f_obs;

//Fill in sample sizes of comps, possibly sampled in nonconsec yrs 
//Used primarily for output in R object   

      nsamp_comm_lenc_allyr=missing;
	  //nsamp_cL_lenc_allyr=missing;   
	  //nsamp_GR_lenc_allyr=missing;	  
	  nsamp_GR_agec_allyr=missing;  
	  
      nfish_comm_lenc_allyr=missing;
	  //nfish_cL_lenc_allyr=missing;  
	  //nfish_GR_lenc_allyr=missing;
	  nfish_GR_agec_allyr=missing;
   
      for (iyear=1; iyear<=nyr_comm_lenc; iyear++)
         {if (nsamp_comm_lenc(iyear)>=minSS_comm_lenc)
           {nsamp_comm_lenc_allyr(yrs_comm_lenc(iyear))=nsamp_comm_lenc(iyear);
            nfish_comm_lenc_allyr(yrs_comm_lenc(iyear))=nfish_comm_lenc(iyear);}}
      //for (iyear=1; iyear<=nyr_cL_lenc; iyear++)
      //   {if (nsamp_cL_lenc(iyear)>=minSS_cL_lenc)
      //     {nsamp_cL_lenc_allyr(yrs_cL_lenc(iyear))=nsamp_cL_lenc(iyear);
      //      nfish_cL_lenc_allyr(yrs_cL_lenc(iyear))=nfish_cL_lenc(iyear);}}			
	  //for (iyear=1; iyear<=nyr_GR_lenc; iyear++)                           
      //   {if (nsamp_GR_lenc(iyear)>=minSS_GR_lenc)
      //      {nsamp_GR_lenc_allyr(yrs_GR_lenc(iyear))=nsamp_GR_lenc(iyear);
      //       nfish_GR_lenc_allyr(yrs_GR_lenc(iyear))=nfish_GR_lenc(iyear);}}
	  for (iyear=1; iyear<=nyr_GR_agec; iyear++)  
         {if (nsamp_GR_agec(iyear)>=minSS_GR_agec)
           {nsamp_GR_agec_allyr(yrs_GR_agec(iyear))=nsamp_GR_agec(iyear);
             nfish_GR_agec_allyr(yrs_GR_agec(iyear))=nfish_GR_agec(iyear);}}  
	  
//fill in Fs for msy and per-recruit analyses
  F_msy(1)=0.0;  
  for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
  for (ff=2;ff<=n_iter_spr;ff++) {F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}


//fill in F's, Catch matrices, and log rec dev with zero's
  F_comm.initialize(); L_comm_num.initialize();
  //F_cL.initialize(); L_cL_num.initialize();
  // F_HB.initialize(); L_HB_num.initialize();
  F_GR.initialize(); L_GR_num.initialize();

  F_comm_out.initialize();
  //F_cL_out.initialize();
  // F_HB_out.initialize();
  F_GR_out.initialize();

  sel_comm.initialize();
  sel_comm_vec.initialize();
  sel_HB.initialize();
  sel_GR.initialize();
  

  //sel_comm_block1.initialize();  
  //sel_comm_block2.initialize();  
  //sel_cL_block1.initialize();    
  //sel_cL_block2.initialize(); 
  sel_HB_block1.initialize();
  sel_HB_block2.initialize();
  sel_GR_block1.initialize();
  sel_GR_block2.initialize();
  
  
  log_rec_dev_output.initialize();  
  log_rec_dev=set_log_rec_dev_vals;
  log_Nage_dev_output.initialize();
  log_Nage_dev=set_log_Nage_dev_vals;
 
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
TOP_OF_MAIN_SECTION
  time(&start);
  arrmblsize=20000000;
  gradient_structure::set_MAX_NVAR_OFFSET(1600);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(2000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(10000);

//>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PROCEDURE_SECTION

///////////////////////////////////////////////////////////////////////////////////////////////////
//MCMC stuff
///////////////////////////////////////////////////////////////////////////////////////////////////
  if (mceval_phase())
   {
//   ofstream mcres("Params.txt",ios::app);
//   ofstream mcres("Fit.txt",ios::app);
//   ofstream mcres("Pop_quant.txt",ios::app);
//   ofstream mcres("NLLikelihood.txt",ios::app);
//   ofstream mcres("Tot_Data_NLL.txt",ios::app);
   ofstream mcres("SDNRs.txt",ios::app);
    if(header==1){

   //Extracting Parameters
   /*
   mcres << "Linf K t0 len_CV_val Linf_L  K_L  t0_L len_cv_val_L Linf_F K_F t0_F len_cv_val_F \
             log_Nage_dev2 log_Nage_dev3 log_Nage_dev4 log_Nage_dev5 log_Nage_dev6 log_Nage_dev7 log_Nage_dev8 log_Nage_dev9 log_Nage_dev10 log_Nage_dev11 log_Nage_dev12 log_Nage_dev13 log_Nage_dev14 log_Nage_dev15 log_Nage_dev16 \
             log_R0 steep rec_sigma R_autocorr \
             log_rec_dev86 log_rec_dev87 log_rec_dev88 log_rec_dev89 log_rec_dev90 log_rec_dev91 log_rec_dev92 log_rec_dev93 log_rec_dev94 log_rec_dev95 log_rec_dev96 log_rec_dev97 log_rec_dev98 log_rec_dev99 log_rec_dev00 log_rec_dev01 log_rec_dev02 log_rec_dev03 log_rec_dev04 log_rec_dev05 log_rec_dev06 log_rec_dev07 log_rec_dev08 log_rec_dev09 log_rec_dev10 log_rec_dev11 log_rec_dev12 log_rec_dev13 log_rec_dev14 log_rec_dev15 log_rec_dev16 \
             log_dm_GR_ac log_dm_comm_lc selpar_A50_comm1 selpar_slope_comm1 selpar_A50_GR1 selpar_slope_GR1 selpar_A50_GR2 selpar_slope_GR2 \
             log_q_HB \
             q_RW_log_dev_HB91 q_RW_log_dev_HB92 q_RW_log_dev_HB93 q_RW_log_dev_HB94 q_RW_log_dev_HB95 q_RW_log_dev_HB96 q_RW_log_dev_HB97 q_RW_log_dev_HB98 q_RW_log_dev_HB99 q_RW_log_dev_HB00 q_RW_log_dev_HB01 q_RW_log_dev_HB02 q_RW_log_dev_HB03 q_RW_log_dev_HB04 q_RW_log_dev_HB05 q_RW_log_dev_HB06 q_RW_log_dev_HB07 q_RW_log_dev_HB08 q_RW_log_dev_HB09 q_RW_log_dev_HB10 q_RW_log_dev_HB11 q_RW_log_dev_HB12 q_RW_log_dev_HB13 q_RW_log_dev_HB14 \
             M_constant \
             log_avg_F_comm \
             log_F_dev_comm86 log_F_dev_comm87 log_F_dev_comm88 log_F_dev_comm89 log_F_dev_comm90 log_F_dev_comm91 log_F_dev_comm92 log_F_dev_comm93 log_F_dev_comm94 log_F_dev_comm95 log_F_dev_comm96 log_F_dev_comm97 log_F_dev_comm98 log_F_dev_comm99 log_F_dev_comm00 log_F_dev_comm01 log_F_dev_comm02 log_F_dev_comm03 log_F_dev_comm04 log_F_dev_comm05 log_F_dev_comm06 log_F_dev_comm07 log_F_dev_comm08 log_F_dev_comm09 log_F_dev_comm10 log_F_dev_comm11 log_F_dev_comm12 log_F_dev_comm13  log_F_dev_comm14 log_F_dev_comm15 log_F_dev_comm16 log_F_dev_comm17 \
             log_avg_F_GR \
             log_F_dev_GR86 log_F_dev_GR87 log_F_dev_GR88 log_F_dev_GR89 log_F_dev_GR90 log_F_dev_GR91 log_F_dev_GR92 log_F_dev_GR93 log_F_dev_GR94 log_F_dev_GR95 log_F_dev_GR96 log_F_dev_GR97 log_F_dev_GR98 log_F_dev_GR99 log_F_dev_GR00 log_F_dev_GR01 log_F_dev_GR02 log_F_dev_GR03 log_F_dev_GR04 log_F_dev_GR05 log_F_dev_GR06 log_F_dev_GR07 log_F_dev_GR08 log_F_dev_GR09 log_F_dev_GR10 log_F_dev_GR11 log_F_dev_GR12 log_F_dev_GR13  log_F_dev_GR14 log_F_dev_GR15 log_F_dev_GR16 log_F_dev_GR17 \
             F_init"<<endl;
    */

   //Data Fit Extraction
   /*
   mcres << "HB_cpue91 HB_cpue92 HB_cpue93 HB_cpue94 HB_cpue95 HB_cpue96 HB_cpue97 HB_cpue98 HB_cpue99 HB_cpue00 HB_cpue01 HB_cpue02 HB_cpue03 HB_cpue04 HB_cpue05 HB_cpue06 HB_cpue07 HB_cpue08 HB_cpue09 HB_cpue10 HB_cpue11 HB_cpue12 HB_cpue13 HB_cpue14 HB_cpue15 \
             CommL_86 CommL_87 CommL_88 CommL_89 CommL_90 CommL_91 CommL_92 CommL_93 CommL_94 CommL_95 CommL_96 CommL_97 CommL_98 CommL_99 CommL_00 CommL_01 CommL_02 CommL_03 CommL_04 CommL_05 CommL_06 CommL_07 CommL_08 CommL_09 CommL_10 CommL_11 CommL_12 CommL_13 CommL_14 CommL_15 CommL_16 CommL_17 \
             GRL_86 GRL_87 GRL_88 GRL_89 GRL_90 GRL_91 GRL_92 GRL_93 GRL_94 GRL_95 GRL_96 GRL_97 GRL_98 GRL_99 GRL_00 GRL_01 GRL_02 GRL_03 GRL_04 GRL_05 GRL_06 GRL_07 GRL_08 GRL_09 GRL_10 GRL_11 GRL_12 GRL_13 GRL_14 GRL_15 GRL_16 GRL_17 \
             lcomp_260 lcomp_290 lcomp_320 lcomp_350 lcomp_380 lcomp_410 lcomp_440 lcomp_470 lcomp_500 lcomp_530 lcomp_560 lcomp_590 lcomp_620 lcomp_650 lcomp_680 lcomp_710 lcomp_740 lcomp_770 lcomp_800 lcomp_830 lcomp_860 \
             lcomp_890 lcomp_920 lcomp_950 lcomp_980 lcomp_1010 lcomp_1040 lcomp_1070 lcomp_1100 lcomp_1130 lcomp_1160 lcomp_1190 lcomp_1220 lcomp_1250 lcomp_1280 lcomp_1310 lcomp_1340 lcomp_1370 lcomp_1400 lcomp_1430 lcomp_1460 \
             lcomp_1490 lcomp_1520 lcomp_1550 lcomp_1580 lcomp_1610 lcomp_1640 lcomp_1670 \
             1986_1 1986_2 1986_3 1986_4 1986_5 1986_6 1986_7 1986_8 1986_9 1986_10 1986_11 1986_12 \
             1987_1 1987_2 1987_3 1987_4 1987_5 1987_6 1987_7 1987_8 1987_9 1987_10 1987_11 1987_12 \
             1989_1 1989_2 1989_3 1989_4 1989_5 1989_6 1989_7 1989_8 1989_9 1989_10 1989_11 1989_12 \
             1990_1 1990_2 1990_3 1990_4 1990_5 1990_6 1990_7 1990_8 1990_9 1990_10 1990_11 1990_12 \
             1991_1 1991_2 1991_3 1991_4 1991_5 1991_6 1991_7 1991_8 1991_9 1991_10 1991_11 1991_12 \
             1992_1 1992_2 1992_3 1992_4 1992_5 1992_6 1992_7 1992_8 1992_9 1992_10 1992_11 1992_12 \
             1995_1 1995_2 1995_3 1995_4 1995_5 1995_6 1995_7 1995_8 1995_9 1995_10 1995_11 1995_12 \
             1996_1 1996_2 1996_3 1996_4 1996_5 1996_6 1996_7 1996_8 1996_9 1996_10 1996_11 1996_12 \
             1997_1 1997_2 1997_3 1997_4 1997_5 1997_6 1997_7 1997_8 1997_9 1997_10 1997_11 1997_12 \
             1999_1 1999_2 1999_3 1999_4 1999_5 1999_6 1999_7 1999_8 1999_9 1999_10 1999_11 1999_12 \
             2000_1 2000_2 2000_3 2000_4 2000_5 2000_6 2000_7 2000_8 2000_9 2000_10 2000_11 2000_12 \
             2001_1 2001_2 2001_3 2001_4 2001_5 2001_6 2001_7 2001_8 2001_9 2001_10 2001_11 2001_12 \
             2002_1 2002_2 2002_3 2002_4 2002_5 2002_6 2002_7 2002_8 2002_9 2002_10 2002_11 2002_12 \
             2005_1 2005_2 2005_3 2005_4 2005_5 2005_6 2005_7 2005_8 2005_9 2005_10 2005_11 2005_12 \
             2006_1 2006_2 2006_3 2006_4 2006_5 2006_6 2006_7 2006_8 2006_9 2006_10 2006_11 2006_12 \
             2007_1 2007_2 2007_3 2007_4 2007_5 2007_6 2007_7 2007_8 2007_9 2007_10 2007_11 2007_12 \
             2008_1 2008_2 2008_3 2008_4 2008_5 2008_6 2008_7 2008_8 2008_9 2008_10 2008_11 2008_12 \
             2009_1 2009_2 2009_3 2009_4 2009_5 2009_6 2009_7 2009_8 2009_9 2009_10 2009_11 2009_12 \
             2010_1 2010_2 2010_3 2010_4 2010_5 2010_6 2010_7 2010_8 2010_9 2010_10 2010_11 2010_12 \
             2011_1 2011_2 2011_3 2011_4 2011_5 2011_6 2011_7 2011_8 2011_9 2011_10 2012_11 2011_12 \
             2012_1 2012_2 2012_3 2012_4 2012_5 2012_6 2012_7 2012_8 2012_9 2012_10 2012_11 2012_12 \
             2013_1 2013_2 2013_3 2013_4 2013_5 2013_6 2013_7 2013_8 2013_9 2013_10 2013_11 2013_12 \
             2014_1 2014_2 2014_3 2014_4 2014_5 2014_6 2014_7 2014_8 2014_9 2014_10 2014_11 2014_12 \
             2015_1 2015_2 2015_3 2015_4 2015_5 2015_6 2015_7 2015_8 2015_9 2015_10 2015_11 2015_12 \
             2016_1 2016_2 2016_3 2016_4 2016_5 2016_6 2016_7 2016_8 2016_9 2016_10 2016_11 2016_12 \
             2017_1 2017_2 2017_3 2017_4 2017_5 2017_6 2017_7 2017_8 2017_9 2017_10 2017_11 2017_12"<<endl;
    */

    //Population Quantities, Not going to bother with labels here
   /*
   mcres << "1986_1 1986_2 1986_3 1986_4 1986_5 1986_6 1986_7 1986_8 1986_9 1986_10 1986_11 1986_12 1986_13 1986_14 1986_15 1986_16 \
             1987_1 1987_2 1987_3 1987_4 1987_5 1987_6 1987_7 1987_8 1987_9 1987_10 1987_11 1987_12 1987_13 1987_14 1987_15 1987_16 \
             1988_1 1988_2 1988_3 1988_4 1988_5 1988_6 1988_7 1988_8 1988_9 1988_10 1988_11 1988_12 1988_13 1988_14 1988_15 1988_16 \
             1989_1 1989_2 1989_3 1989_4 1989_5 1989_6 1989_7 1989_8 1989_9 1989_10 1989_11 1989_12 1989_13 1989_14 1989_15 1989_16 \
             1990_1 1990_2 1990_3 1990_4 1990_5 1990_6 1990_7 1990_8 1990_9 1990_10 1990_11 1990_12 1990_13 1990_14 1990_15 1990_16 \
             1991_1 1991_2 1991_3 1991_4 1991_5 1991_6 1991_7 1991_8 1991_9 1991_10 1991_11 1991_12 1991_13 1991_14 1991_15 1991_16 \
             1992_1 1992_2 1992_3 1992_4 1992_5 1992_6 1992_7 1992_8 1992_9 1992_10 1992_11 1992_12 1992_13 1992_14 1992_15 1992_16 \
             1993_1 1993_2 1993_3 1993_4 1993_5 1993_6 1993_7 1993_8 1993_9 1993_10 1993_11 1993_12 1993_13 1993_14 1993_15 1993_16 \
             1994_1 1994_2 1994_3 1994_4 1994_5 1994_6 1994_7 1994_8 1994_9 1994_10 1994_11 1994_12 1994_13 1994_14 1994_15 1994_16 \
             1995_1 1995_2 1995_3 1995_4 1995_5 1995_6 1995_7 1995_8 1995_9 1995_10 1995_11 1995_12 1995_13 1995_14 1995_15 1995_16 \
             1996_1 1996_2 1996_3 1996_4 1996_5 1996_6 1996_7 1996_8 1996_9 1996_10 1996_11 1996_12 1996_13 1996_14 1996_15 1996_16 \
             1997_1 1997_2 1997_3 1997_4 1997_5 1997_6 1997_7 1997_8 1997_9 1997_10 1997_11 1997_12 1997_13 1997_14 1997_15 1997_16 \
             1998_1 1998_2 1998_3 1998_4 1998_5 1998_6 1998_7 1998_8 1998_9 1998_10 1998_11 1998_12 1998_13 1998_14 1998_15 1998_16 \
             1999_1 1999_2 1999_3 1999_4 1999_5 1999_6 1999_7 1999_8 1999_9 1999_10 1999_11 1999_12 1999_13 1999_14 1999_15 1999_16 \
             2000_1 2000_2 2000_3 2000_4 2000_5 2000_6 2000_7 2000_8 2000_9 2000_10 2000_11 2000_12 2000_13 2000_14 2000_15 2000_16 \
             2001_1 2001_2 2001_3 2001_4 2001_5 2001_6 2001_7 2001_8 2001_9 2001_10 2001_11 2001_12 2001_13 2001_14 2001_15 2001_16 \
             2002_1 2002_2 2002_3 2002_4 2002_5 2002_6 2002_7 2002_8 2002_9 2002_10 2002_11 2002_12 2002_13 2002_14 2002_15 2002_16 \
             2003_1 2003_2 2003_3 2003_4 2003_5 2003_6 2003_7 2003_8 2003_9 2003_10 2003_11 2003_12 2003_13 2003_14 2003_15 2003_16 \
             2004_1 2004_2 2004_3 2004_4 2004_5 2004_6 2004_7 2004_8 2004_9 2004_10 2004_11 2004_12 2004_13 2004_14 2004_15 2004_16 \
             2005_1 2005_2 2005_3 2005_4 2005_5 2005_6 2005_7 2005_8 2005_9 2005_10 2005_11 2005_12 2005_13 2005_14 2005_15 2005_16 \
             2006_1 2006_2 2006_3 2006_4 2006_5 2006_6 2006_7 2006_8 2006_9 2006_10 2006_11 2006_12 2006_13 2006_14 2006_15 2006_16 \
             2007_1 2007_2 2007_3 2007_4 2007_5 2007_6 2007_7 2007_8 2007_9 2007_10 2007_11 2007_12 2007_13 2007_14 2007_15 2007_16 \
             2008_1 2008_2 2008_3 2008_4 2008_5 2008_6 2008_7 2008_8 2008_9 2008_10 2008_11 2008_12 2008_13 2008_14 2008_15 2008_16 \
             2009_1 2009_2 2009_3 2009_4 2009_5 2009_6 2009_7 2009_8 2009_9 2009_10 2009_11 2009_12 2009_13 2009_14 2009_15 2009_16 \
             2010_1 2010_2 2010_3 2010_4 2010_5 2010_6 2010_7 2010_8 2010_9 2010_10 2010_11 2010_12 2010_13 2010_14 2010_15 2010_16 \
             2011_1 2011_2 2011_3 2011_4 2011_5 2011_6 2011_7 2011_8 2011_9 2011_10 2012_11 2011_12 2011_13 2011_14 2011_15 2011_16 \
             2012_1 2012_2 2012_3 2012_4 2012_5 2012_6 2012_7 2012_8 2012_9 2012_10 2012_11 2012_12 2012_13 2012_14 2012_15 2012_16 \
             2013_1 2013_2 2013_3 2013_4 2013_5 2013_6 2013_7 2013_8 2013_9 2013_10 2013_11 2013_12 2013_13 2013_14 2013_15 2013_16 \
             2014_1 2014_2 2014_3 2014_4 2014_5 2014_6 2014_7 2014_8 2014_9 2014_10 2014_11 2014_12 2014_13 2014_14 2014_15 2014_16 \
             2015_1 2015_2 2015_3 2015_4 2015_5 2015_6 2015_7 2015_8 2015_9 2015_10 2015_11 2015_12 2015_13 2015_14 2015_15 2015_16 \
             2016_1 2016_2 2016_3 2016_4 2016_5 2016_6 2016_7 2016_8 2016_9 2016_10 2016_11 2016_12 2016_13 2016_14 2016_15 2016_16 \
             2017_1 2017_2 2017_3 2017_4 2017_5 2017_6 2017_7 2017_8 2017_9 2017_10 2017_11 2017_12 2017_13 2017_14 2017_15 2017_16 \
             2018_1 2018_2 2018_3 2018_4 2018_5 2018_6 2018_7 2018_8 2018_9 2018_10 2018_11 2018_12 2018_13 2018_14 2018_15 2018_16 \
             Wt_1 Wt_2 Wt_3 Wt_4 Wt_5 Wt_6 Wt_7 Wt_8 Wt_9 Wt_10 Wt_11 Wt_12 Wt_13 Wt_14 Wt_15 Wt_16 \
             SB_1986 SB_1987 SB_1988 SB_1989 SB_1990 SB_1991 SB_1992 SB_1993 SB_1994 SB_1995 SB_1996 SB_1997 SB_1998 SB_1999 SB_2000 SB_2001 SB_2002 SB_2003 SB_2004 SB_2005 SB_2006 SB_2007 SB_2008 SB_2009 SB_2010 SB_2011 SB_2012 SB_2013 SB_2014 SB_2015 SB_2016 SB_2017 \
             Llb_1986 Llb_1987 Llb_1988 Llb_1989 Llb_1990 Llb_1991 Llb_1992 Llb_1993 Llb_1994 Llb_1995 Llb_1996 Llb_1997 Llb_1998 Llb_1999 Llb_2000 Llb_2001 Llb_2002 Llb_2003 Llb_2004 Llb_2005 Llb_2006 Llb_2007 Llb_2008 Llb_2009 Llb_2010 Llb_2011 Llb_2012 Llb_2013 Llb_2014 Llb_2015 Llb_2016 Llb_2017 \
             L_1986 L_1987 L_1988 L_1989 L_1990 L_1991 L_1992 L_1993 L_1994 L_1995 L_1996 L_1997 L_1998 L_1999 L_2000 L_2001 L_2002 L_2003 L_2004 L_2005 L_2006 L_2007 L_2008 L_2009 L_2010 L_2011 L_2012 L_2013 L_2014 L_2015 L_2016 L_2017"<<endl;
    */
    header=0;
    }

   //Parameter Extraction
   /*
   mcres << Linf << " " << K << " " << t0 << " " << len_cv_val << " " << Linf_L << " " << K_L << " " << t0_L << " " << len_cv_val_L << " " << Linf_F << " " << K_F << " " << t0_F << " " << len_cv_val_F << " "
         << log_Nage_dev << " " << log_R0 << " " << steep << " " << rec_sigma << " " << R_autocorr << " " << log_rec_dev << " "
         << log_dm_GR_ac << " " << log_dm_comm_lc << " " << selpar_A50_comm1 << " " << selpar_slope_comm1 << " " << selpar_A50_GR1 << " " << selpar_slope_GR1 << " " << selpar_A50_GR2 << " " << selpar_slope_GR2 << " "
         << log_q_HB << " " << q_RW_log_dev_HB << " " << M_constant << " " << log_avg_F_comm << " " << log_F_dev_comm << " " << log_avg_F_GR << " " << log_F_dev_GR << " " << F_init <<  endl;
   */

   //Data Fit Extraction
   /*
   mcres <<pred_HB_cpue<<" "
         <<pred_comm_L_klb(styr_comm_L,endyr_comm_L)<<" "
         <<pred_GR_L_knum(styr_GR_L,endyr_GR_L)<<" "
         <<pred_comm_lenc<<" "
         <<pred_GR_agec(1)<<" "<<pred_GR_agec(2)<<" "<<pred_GR_agec(3)<<" "<<pred_GR_agec(4)<<" "<<pred_GR_agec(5)<<" "<<pred_GR_agec(6)<<" "<<pred_GR_agec(7)<<" "<<pred_GR_agec(8)<<" "<<pred_GR_agec(9)<<" "<<pred_GR_agec(10)<<" "
         <<pred_GR_agec(11)<<" "<<pred_GR_agec(12)<<" "<<pred_GR_agec(13)<<" "<<pred_GR_agec(14)<<" "<<pred_GR_agec(15)<<" "<<pred_GR_agec(16)<<" "<<pred_GR_agec(17)<<" "<<pred_GR_agec(18)<<" "<<pred_GR_agec(19)<<" "<<pred_GR_agec(20)<<" "
         <<pred_GR_agec(21)<<" "<<pred_GR_agec(22)<<" "<<pred_GR_agec(23)<<" "<<pred_GR_agec(24)<<" "<<pred_GR_agec(25)<<" "<<pred_GR_agec(26)<<endl;
   */

   //Population Quantities,
   /*
   mcres <<N(1986)<<" "<<N(1987)<<" "<<N(1988)<<" "<<N(1989)<<" "<<N(1990)<<" "<<N(1991)<<" "<<N(1992)<<" "<<N(1993)<<" "<<N(1994)<<" "<<N(1995)<<" "<<N(1996)<<" "<<N(1997)<<" "<<N(1998)<<" "<<N(1999)<<" "<<N(2000)<<" "<<N(2001)<<" "<<N(2002)<<" "<<N(2003)<<" "<<N(2004)<<" "<<N(2005)<<" "<<N(2006)<<" "<<N(2007)<<" "<<N(2008)<<" "<<N(2009)<<" "<<N(2010)<<" "<<N(2011)<<" "<<N(2012)<<" "<<N(2013)<<" "<<N(2014)<<" "<<N(2015)<<" "<<N(2016)<<" "<<N(2017)<<" "<<N(2018)<<" "
         <<wgt_lb<<" "
         <<SSB<<" "
         <<pred_comm_L_klb(styr,endyr)+pred_GR_L_klb(styr,endyr)<<" "
         <<pred_comm_L_knum(styr,endyr)+pred_GR_L_knum(styr,endyr)<<endl;
   */
   
   //Negative Log Likelihood of each data point
   //mcres << LL_HB_cpue << " " << LL_comm_L << " " << LL_GR_L << " " << LL_comm_lc << " " << LL_GR_ac << endl;
   //Tot NLL (as a check to the above)
   //mcres <<f_HB_cpue<<" "<<f_comm_L<<" "<<f_GR_L<<" "<<f_comm_lenc<<" "<<f_GR_agec<<endl;
   mcres<<sdnr_lognormal(pred_HB_cpue, obs_HB_cpue, HB_cpue_cv, w_I_HB)<<endl;
   mcres.close();
   }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //cout<<"start"<<endl;
 
 //get_M_at_age(); //Needed only if M is estimated
 
  get_length_weight_at_age(); 
  //cout << "got length, weight, fecundity transitions" <<endl;
  get_reprod();
  //cout << "got reprod" << endl;
  get_length_at_age_dist(); 
  //cout<< "got predicted length at age distribution"<<endl;
  get_weight_at_age_landings();
  //cout<< "got weight at age of landings"<<endl; 
  get_spr_F0();
  //cout << "got F0 spr" << endl;
  get_selectivity(); 
  //cout << "got selectivity" << endl;
  get_mortality(); 
  // cout << "got mortalities" << endl;
  get_bias_corr(); 
  //cout<< "got recruitment bias correction" << endl;
  get_numbers_at_age(); 
  //cout << "got numbers at age" << endl;
  get_landings_numbers();
  //cout << "got landings in numbers" << endl;
  get_landings_wgt();
  //cout << "got landings in wgt" << endl;
  // get_dead_discards(); 
  //cout << "got dead discards in num and wgt" << endl;
  get_catchability_fcns(); 
  //cout << "got catchability_fcns" << endl;
  get_indices();
  //cout << "got indices" << endl;
  get_length_comps();
  // cout<< "got length comps"<< endl;
   get_age_comps();
  //cout<< "got age comps"<< endl;
   evaluate_objective_function();
   //cout << "objective function calculations complete" << endl;
 
 
FUNCTION get_length_weight_at_age
	//population total length in mm
    //compute mean length (mm TL) and weight (whole) at age
    meanlen_TL=Linf*(1.0-mfexp(-K*(agebins-t0+0.5))); //Actually fork length     
    wgt_kg=wgtpar_a*pow(meanlen_TL,wgtpar_b);             //whole wgt in kg 
    wgt_g=wgt_kg/g2kg;                                    //convert wgt in kg to weight in g    
    wgt_mt=wgt_g*g2mt;                                    //convert weight in g to weight in mt
    wgt_klb=mt2klb*wgt_mt;                                //1000 lb of whole wgt
    wgt_lb=mt2lb*wgt_mt;                                  //lb of whole wgt
    
	//All fisheries
	meanlen_TL_L=Linf_L*(1.0-mfexp(-K_L*(agebins-t0_L+0.5)));    //Landings total length in mm
    wgt_kg_L=wgtpar_a*pow(meanlen_TL_L,wgtpar_b);             //whole wgt in kg 
    wgt_g_L=wgt_kg_L/g2kg;                                    //convert wgt in kg to weight in g    
    wgt_mt_L=wgt_g_L*g2mt;                                    //convert weight in g to weight in mt
    wgt_klb_L=mt2klb*wgt_mt_L;                               //1000 lb of whole wgt
    wgt_lb_L=mt2lb*wgt_mt_L;                                 //1000 lb of whole wgt
	
	//Females
	meanlen_TL_F=Linf_F*(1.0-mfexp(-K_F*(agebins-t0_F+0.5)));    //Landings total length in mm
    wgt_kg_F=wgtpar_a*pow(meanlen_TL_F,wgtpar_b);             //whole wgt in kg 
    wgt_g_F=wgt_kg_F/g2kg;                                    //convert wgt in kg to weight in g    
    wgt_mt_F=wgt_g_F*g2mt;                                    //convert weight in g to weight in mt
    wgt_klb_F=mt2klb*wgt_mt_F;                               //1000 lb of whole wgt
    wgt_lb_F=mt2lb*wgt_mt_F;                                 //1000 lb of whole wgt
	
	//batchfec = mfexp(batchfecpar_a + batchfecpar_b*meanlen_TL);  // batch fecundity at length [should be batchfec = exp(a+bL) based on Harris 2004]
	//fec = batchfec*nbatch/fecpar_scale;                          // annual fecundity at length scaled to fecpar_scale units
	
FUNCTION get_reprod 
	
    //reprod=elem_prod(prop_f,elem_prod(maturity_f,fec));
	reprod=elem_prod(elem_prod(prop_f,maturity_f),wgt_mt_F);
	reprodknum=elem_prod(prop_f,maturity_f)/1000.0;
	//+elem_prod(prop_m,maturity_m)),wg_mt);
 
FUNCTION get_length_at_age_dist
  //compute matrix of length at age, based on the normal distribution
    //population
	for (iage=1;iage<=nages;iage++)
   {len_cv(iage)=len_cv_val;
    len_sd(iage)=meanlen_TL(iage)*len_cv(iage);
	zscore_lzero=(0.0-meanlen_TL(iage))/len_sd(iage); 
	cprob_lzero=cumd_norm(zscore_lzero);
	
    //All fishery dependent
	//len_cv_L(iage)=mfexp(log_len_cv_L+log_len_cv_dev_L(iage));
    len_cv_L(iage)=len_cv_val_L;
    len_sd_L(iage)=meanlen_TL_L(iage)*len_cv_L(iage);
	zscore_lzero_L=(0.0-meanlen_TL_L(iage))/len_sd_L(iage); 
	cprob_lzero_L=cumd_norm(zscore_lzero_L);    
	
	//Females
	//len_cv_L(iage)=mfexp(log_len_cv_L+log_len_cv_dev_L(iage));
    len_cv_F(iage)=len_cv_val_F;
    len_sd_F(iage)=meanlen_TL_F(iage)*len_cv_F(iage);
	zscore_lzero_F=(0.0-meanlen_TL_F(iage))/len_sd_F(iage); 
	cprob_lzero_F=cumd_norm(zscore_lzero_F);
	
    //first length bin
	//population
    zscore_len=((lenbins(1)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage);
    cprob_lenvec(1)=cumd_norm(zscore_len);          //includes any probability mass below zero
    lenprob(iage,1)=cprob_lenvec(1)-cprob_lzero;    //removes any probability mass below zero
	
	//All fishery dependent
    zscore_len_L=((lenbins(1)+0.5*lenbins_width)-meanlen_TL_L(iage)) / len_sd_L(iage);
    cprob_lenvec_L(1)=cumd_norm(zscore_len_L);          //includes any probability mass below zero
    lenprob_L(iage,1)=cprob_lenvec_L(1)-cprob_lzero_L;    //removes any probability mass below zero
	
	//Females
    zscore_len_F=((lenbins(1)+0.5*lenbins_width)-meanlen_TL_F(iage)) / len_sd_F(iage);
    cprob_lenvec_F(1)=cumd_norm(zscore_len_F);          //includes any probability mass below zero
    lenprob_F(iage,1)=cprob_lenvec_F(1)-cprob_lzero_F;    //removes any probability mass below zero
	
    //most other length bins  
    //population
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len=((lenbins(ilen)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage); 
		cprob_lenvec(ilen)=cumd_norm(zscore_len);
        lenprob(iage,ilen)=cprob_lenvec(ilen)-cprob_lenvec(ilen-1);
      }
	
	//All fishery dependent	
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len_L=((lenbins(ilen)+0.5*lenbins_width)-meanlen_TL_L(iage)) / len_sd_L(iage); 
		cprob_lenvec_L(ilen)=cumd_norm(zscore_len_L);
        lenprob_L(iage,ilen)=cprob_lenvec_L(ilen)-cprob_lenvec_L(ilen-1);
      }
	  
	//Females	
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len_F=((lenbins(ilen)+0.5*lenbins_width)-meanlen_TL_F(iage)) / len_sd_F(iage); 
		cprob_lenvec_F(ilen)=cumd_norm(zscore_len_F);
        lenprob_F(iage,ilen)=cprob_lenvec_F(ilen)-cprob_lenvec_F(ilen-1);
      }
	  
    //last length bin is a plus group
	//population
    zscore_len=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage); 
	lenprob(iage,nlenbins)=1.0-cumd_norm(zscore_len);
      lenprob(iage)=lenprob(iage)/(1.0-cprob_lzero);  //renormalize to account for any prob mass below size=0
	  
	//All fishery dependent
    zscore_len_L=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_TL_L(iage)) / len_sd_L(iage); 
	lenprob_L(iage,nlenbins)=1.0-cumd_norm(zscore_len_L);
      lenprob_L(iage)=lenprob_L(iage)/(1.0-cprob_lzero_L);  //renormalize to account for any prob mass below size=0

	//Females
    zscore_len_F=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_TL_F(iage)) / len_sd_F(iage); 
	lenprob_F(iage,nlenbins)=1.0-cumd_norm(zscore_len_F);
      lenprob_F(iage)=lenprob_F(iage)/(1.0-cprob_lzero_F);  //renormalize to account for any prob mass below size=0
   }
   
  //fleet and survey specific length probs, all assumed here to equal the popn
  lenprob_comm=lenprob_L;
  //lenprob_cL=lenprob; 
  lenprob_HB=lenprob;
  //lenprob_GR=lenprob;  
  
FUNCTION get_weight_at_age_landings  ///***in whole weight
  
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    len_comm_mm(iyear)=meanlen_TL_L;  
    wholewgt_comm_klb(iyear)=wgt_klb_L; 
    //len_cL_mm(iyear)=meanlen_TL;  
    //wholewgt_cL_klb(iyear)=wgt_klb; 	
    len_HB_mm(iyear)=meanlen_TL_L;
    wholewgt_HB_klb(iyear)=wgt_klb_L;
    len_GR_mm(iyear)=meanlen_TL_L;
    wholewgt_GR_klb(iyear)=wgt_klb_L;
  }  
 
 
FUNCTION get_spr_F0
  //at mdyr, apply half this yr's mortality, half next yr's
  N_spr_F0(1)=1.0*mfexp(-1.0*M(1)*spawn_time_frac); //at peak spawning time
  N_bpr_F0(1)=1.0;      //at start of year
  for (iage=2; iage<=nages; iage++)
  { N_spr_F0(iage)=N_spr_F0(iage-1)*mfexp(-1.0*(M(iage-1)*(1.0-spawn_time_frac) + M(iage)*spawn_time_frac)); 
    N_bpr_F0(iage)=N_bpr_F0(iage-1)*mfexp(-1.0*(M(iage-1)));    
  }
  N_spr_F0(nages)=N_spr_F0(nages)/(1.0-mfexp(-1.0*M(nages))); //plus group (sum of geometric series)
  N_bpr_F0(nages)=N_bpr_F0(nages)/(1.0-mfexp(-1.0*M(nages)));
  
  spr_F0=sum(elem_prod(N_spr_F0,reprod)); 
  bpr_F0=sum(elem_prod(N_bpr_F0,wgt_mt));    



FUNCTION get_selectivity
  sel_comm_vec=logistic(agebins, selpar_A50_comm1, selpar_slope_comm1);
  sel_GR_block1=logistic(agebins, selpar_A50_GR1, selpar_slope_GR1);
  sel_GR_block2=logistic(agebins, selpar_A50_GR2, selpar_slope_GR2);
  sel_HB_block1=sel_GR_block1; // Use GR selectivity for HB
  sel_HB_block2=sel_GR_block1; // Use GR selectivity for HB
   
  //-------- comm --------// 
  for (iyear=styr; iyear<=endyr; iyear++)
	{sel_comm(iyear) = sel_comm_vec;}
  
  //---- GR and HB ----//   
  //BLOCK 1 for selex  
  for (iyear=styr; iyear<=endyr_selex_phase1_GR; iyear++)
   {     
    sel_HB(iyear)=sel_HB_block1;
    sel_GR(iyear)=sel_GR_block1;
   }
  //BLOCK 2 for selex  
  for (iyear=(endyr_selex_phase1_GR+1); iyear<=endyr; iyear++)//iyear<=endyr_selex_phase2_GR; iyear++)
   {    
    sel_HB(iyear)=sel_HB_block2;
    sel_GR(iyear)=sel_GR_block2;
   }
   
   
FUNCTION get_mortality
  Fsum.initialize();
  Fapex.initialize();
  F.initialize();
  //initialization F is avg from first 3 yrs of observed landings
  log_F_dev_init_comm=sum(log_F_dev_comm(styr_comm_L,(styr_comm_L+2)))/3.0;  
  //log_F_dev_init_cL=sum(log_F_dev_cL(styr_cL_L,(styr_cL_L+2)))/3.0;          
  log_F_dev_init_GR=sum(log_F_dev_GR(styr_GR_L,(styr_GR_L+2)))/3.0;         
  
  for (iyear=styr; iyear<=endyr; iyear++) 
  {
    if(iyear>=styr_comm_L & iyear<=endyr_comm_L) //spans full time series
		{F_comm_out(iyear)=mfexp(log_avg_F_comm+log_F_dev_comm(iyear));}     
    F_comm(iyear)=sel_comm(iyear)*F_comm_out(iyear);
    Fsum(iyear)+=F_comm_out(iyear);
    

	if(iyear>=styr_GR_L & iyear<=endyr_GR_L) //starts in 1981
		{F_GR_out(iyear)=mfexp(log_avg_F_GR+log_F_dev_GR(iyear));}    
    if (iyear<styr_GR_L)
		{F_GR_out(iyear)=mfexp(log_avg_F_GR+log_F_dev_init_GR);}
	F_GR(iyear)=sel_GR(iyear)*F_GR_out(iyear); 
    Fsum(iyear)+=F_GR_out(iyear);
 
    //Total F at age
    F(iyear)=F_comm(iyear);  //first in additive series (NO +=)
    //F(iyear)+=F_cL(iyear);
    // F(iyear)+=F_HB(iyear);
    F(iyear)+=F_GR(iyear);
    
    Fapex(iyear)=max(F(iyear));
    Z(iyear)=M+F(iyear);
    
   }  //end iyear 

FUNCTION get_bias_corr
  var_rec_dev=norm2(log_rec_dev(styr_rec_dev,endyr_rec_dev)-
              sum(log_rec_dev(styr_rec_dev,endyr_rec_dev))/nyrs_rec)
              /(nyrs_rec-1.0);                           
  //if (set_BiasCor <= 0.0) {BiasCor=mfexp(var_rec_dev/2.0);}   //bias correction based on empirical residuals
  rec_sigma_sq=square(rec_sigma);
  if (set_BiasCor <= 0.0) {BiasCor=mfexp(rec_sigma_sq/2.0);}   //bias correction based on Rsigma               
  else {BiasCor=set_BiasCor;}

FUNCTION get_numbers_at_age
//Initialization
  R0=mfexp(log_R0);
  S0=spr_F0*R0;
  //R_virgin=SR_eq_func(R0, steep, spr_F0, spr_F0, BiasCor, SR_switch);
  R_virgin=BiasCor*R0; //changed to move away from an SR relationship
  B0=bpr_F0*R_virgin;   
  B0_q_DD=R_virgin*sum(elem_prod(N_bpr_F0(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages))); 
  
  // Commented out code block from Erik's ASPM for red snapper  
  F_init_denom=mfexp(log_avg_F_comm+log_F_dev_init_comm)+mfexp(log_avg_F_GR+log_F_dev_init_GR); //+mfexp(log_avg_F_cL+log_F_dev_init_cL)
  F_init_comm_prop= mfexp(log_avg_F_comm+log_F_dev_init_comm)/F_init_denom;
  //F_init_cL_prop= mfexp(log_avg_F_cL+log_F_dev_init_cL)/F_init_denom;
  F_init_GR_prop= mfexp(log_avg_F_GR+log_F_dev_init_GR)/F_init_denom;
  
  F_initial=sel_comm(styr)*F_init*F_init_comm_prop+
            //sel_cL(styr)*F_init*F_init_cL_prop+
            sel_GR(styr)*F_init*F_init_GR_prop;
			
	//F_initial=sel_initial*F_init;		
  Z_initial=M+F_initial;  

//Initial equilibrium age structure
  N_spr_initial(1)=1.0*mfexp(-1.0*Z_initial(1)*spawn_time_frac); //at peak spawning time;
  for (iage=2; iage<=nages; iage++)
    {
      N_spr_initial(iage)=N_spr_initial(iage-1)*
                   mfexp(-1.0*(Z_initial(iage-1)*(1.0-spawn_time_frac) + Z_initial(iage)*spawn_time_frac)); 
    }
  N_spr_initial(nages)=N_spr_initial(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  spr_initial=sum(elem_prod(N_spr_initial,reprod));
  //if (styr==styr_rec_dev) {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, 1.0, SR_switch);} //without bias correction (deviation added later)
  //else {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, BiasCor, SR_switch);} //with bias correction
  if (styr==styr_rec_dev) {R1=R0;} //without bias correction (deviation added later)
  else {R1=BiasCor*R0;} //with bias correction

  if(R1<10.0) {R1=10.0;} //Avoid unrealistically low popn sizes during search algorithm

//Compute equilibrium age structure for first year
  N_initial_eq(1)=R1;
  for (iage=2; iage<=nages; iage++)
  {
    N_initial_eq(iage)=N_initial_eq(iage-1)*
        mfexp(-1.0*(Z_initial(iage-1)));    
  }
  //plus group calculation
  N_initial_eq(nages)=N_initial_eq(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  
//Add deviations to initial equilibrium N
  N(styr)(2,nages)=elem_prod(N_initial_eq(2,nages),mfexp(log_Nage_dev));
   
  if (styr==styr_rec_dev) {N(styr,1)=N_initial_eq(1)*mfexp(log_rec_dev(styr_rec_dev));}
  else {N(styr,1)=N_initial_eq(1);}
  
  N_mdyr(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*0.5))); //mid year 
  N_spawn(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*spawn_time_frac))); //peak spawning time 

  SSB(styr)=sum(elem_prod(N_spawn(styr),reprod));
  SSB_knum(styr)=sum(elem_prod(N_spawn(styr),reprodknum));
  B_q_DD(styr)=sum(elem_prod(N(styr)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
    
//Rest of years 
  for (iyear=styr; iyear<endyr; iyear++)
  {
    if(iyear<(styr_rec_dev-1)||iyear>(endyr_rec_dev-1)) //recruitment follows S-R curve (with bias correction) exactly
    {
		N(iyear+1,1)=BiasCor*R0; //Changed to use ave rec instead of SR relationship
		//N(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch);
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));
		SSB_knum(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprodknum));
		B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));       
    }
    else   //recruitment follows S-R curve with lognormal deviation
    {
        N(iyear+1,1)=R0*mfexp(log_rec_dev(iyear+1)); //Changed to use ave rec instead of SR relationship
        //N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch)*mfexp(log_rec_dev(iyear+1));
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));
		SSB_knum(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprodknum));
        B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
    }
  }
  
  //last year (projection) has no recruitment variability
  N(endyr+1,1)=BiasCor*R0; //Changed to use ave rec instead of SR relationship
  //N(endyr+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(endyr),SR_switch);
  N(endyr+1)(2,nages)=++elem_prod(N(endyr)(1,nages-1),(mfexp(-1.*Z(endyr)(1,nages-1))));
  N(endyr+1,nages)+=N(endyr,nages)*mfexp(-1.*Z(endyr,nages)); //plus group

  
FUNCTION get_landings_numbers //Baranov catch eqn
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      L_comm_num(iyear,iage)=N(iyear,iage)*F_comm(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      //L_cL_num(iyear,iage)=N(iyear,iage)*F_cL(iyear,iage)*
        //(1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);		
      L_GR_num(iyear,iage)=N(iyear,iage)*F_GR(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);        
    }          
    pred_comm_L_knum(iyear)=sum(L_comm_num(iyear))/1000.0;
	//pred_cL_L_knum(iyear)=sum(L_cL_num(iyear))/1000.0;
    pred_GR_L_knum(iyear)=sum(L_GR_num(iyear))/1000.0;
  }

 
FUNCTION get_landings_wgt
  for (iyear=styr; iyear<=endyr; iyear++)
  {    
    L_comm_klb(iyear)=elem_prod(L_comm_num(iyear),wholewgt_comm_klb(iyear));     //in 1000 lb whole weight
	//L_cL_klb(iyear)=elem_prod(L_cL_num(iyear),wholewgt_cL_klb(iyear));     //in 1000 lb whole weight
    // L_HB_klb(iyear)=elem_prod(L_HB_num(iyear),wholewgt_HB_klb(iyear));     //in 1000 lb whole weight
    L_GR_klb(iyear)=elem_prod(L_GR_num(iyear),wholewgt_GR_klb(iyear));     //in 1000 lb whole weight
    
    pred_comm_L_klb(iyear)=sum(L_comm_klb(iyear));
	//pred_cL_L_klb(iyear)=sum(L_cL_klb(iyear));
    // pred_HB_L_klb(iyear)=sum(L_HB_klb(iyear));
    pred_GR_L_klb(iyear)=sum(L_GR_klb(iyear));    
  }
   
FUNCTION get_catchability_fcns    
 //Get rate increase if estimated, otherwise fixed above
  if (set_q_rate_phase>0.0)
  {
  	  
      for (iyear=styr_HB_cpue; iyear<=endyr_HB_cpue; iyear++)
      {   if (iyear>styr_HB_cpue & iyear <=2003) 
          {//q_rate_fcn_HB(iyear)=(1.0+q_rate)*q_rate_fcn_HB(iyear-1); //compound
             q_rate_fcn_HB(iyear)=(1.0+(iyear-styr_HB_cpue)*q_rate)*q_rate_fcn_HB(styr_HB_cpue);  //linear
          }
          if (iyear>2003) {q_rate_fcn_HB(iyear)=q_rate_fcn_HB(iyear-1);} 
      }      
  } //end q_rate conditional      

 //Get density dependence scalar (=1.0 if density independent model is used)   
  if (q_DD_beta>0.0) 
  {
    B_q_DD+=dzero;
    for (iyear=styr;iyear<=endyr;iyear++)
        {q_DD_fcn(iyear)=pow(B0_q_DD,q_DD_beta)*pow(B_q_DD(iyear),-q_DD_beta);}
          //{q_DD_fcn(iyear)=1.0+4.0/(1.0+mfexp(0.75*(B_q_DD(iyear)-0.1*B0_q_DD))); }
  }  

     
FUNCTION get_indices
//---Predicted CPUEs------------------------

 //HB  cpue
  q_HB(styr_HB_cpue)=mfexp(log_q_HB); 
  for (iyear=styr_HB_cpue; iyear<=endyr_HB_cpue; iyear++)
  {   
      N_HB(iyear)=elem_prod(N_mdyr(iyear),sel_HB(iyear)); 
      pred_HB_cpue(iyear)=q_HB(iyear)*q_rate_fcn_HB(iyear)*q_DD_fcn(iyear)*sum(N_HB(iyear));
      if (iyear<endyr_HB_cpue){q_HB(iyear+1)=q_HB(iyear)*mfexp(q_RW_log_dev_HB(iyear));}
  }
  
FUNCTION get_length_comps

  //comm lines 
  
  for (iyear=1;iyear<=nyr_comm_lenc_pool;iyear++)
   {pred_comm_lenc_yr(iyear)=(L_comm_num(yrs_comm_lenc_pool(iyear))*lenprob_comm)/sum(L_comm_num(yrs_comm_lenc_pool(iyear)));}
  
  pred_comm_lenc.initialize();
  for (iyear=1;iyear<=nyr_comm_lenc_pool;iyear++)
    {pred_comm_lenc(1) += nfish_comm_lenc_pool(iyear) * pred_comm_lenc_yr(iyear);}
  pred_comm_lenc(1)=pred_comm_lenc(1)/sum(nfish_comm_lenc_pool);
  
  ////comm longline 
  //for (iyear=1;iyear<=nyr_cL_lenc;iyear++)
  //{pred_cL_lenc(iyear)=(L_cL_num(yrs_cL_lenc(iyear))*lenprob_cL)/sum(L_cL_num(yrs_cL_lenc(iyear)));}

  //general rec   
  //for (iyear=1;iyear<=nyr_GR_lenc;iyear++) 
  //{pred_GR_lenc(iyear)=(L_GR_num(yrs_GR_lenc(iyear))*lenprob_GR)/sum(L_GR_num(yrs_GR_lenc(iyear)));}

FUNCTION get_age_comps

  //Recreational 
  for (iyear=1;iyear<=nyr_GR_agec;iyear++)
  {
    ErrorFree_GR_agec(iyear)=L_GR_num(yrs_GR_agec(iyear))/sum(L_GR_num(yrs_GR_agec(iyear)));
    pred_GR_agec_allages(iyear)=age_error*ErrorFree_GR_agec(iyear); 
    for (iage=1; iage<=nages_agec; iage++) {pred_GR_agec(iyear,iage)=pred_GR_agec_allages(iyear,iage);} 
    //for (iage=(nages_agec+1); iage<=nages; iage++) {pred_GR_agec(iyear,nages_agec)+=pred_GR_agec_allages(iyear,iage);} //plus group                        
  }
  
////--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_weighted_current 
  F_temp_sum=0.0;
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_comm+
        sum(log_F_dev_comm((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
    
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_GR+
        sum(log_F_dev_GR((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);

  F_comm_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_comm+
        sum(log_F_dev_comm((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  
  F_GR_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_GR+
        sum(log_F_dev_GR((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  
  log_F_dev_end_comm=sum(log_F_dev_comm((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  //log_F_dev_end_cL=sum(log_F_dev_cL((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted; 
  log_F_dev_end_GR=sum(log_F_dev_GR((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  

  F_end_L=sel_comm(endyr)*mfexp(log_avg_F_comm+log_F_dev_end_comm)+
          //sel_cL(endyr)*mfexp(log_avg_F_cL+log_F_dev_end_cL)+
		  sel_GR(endyr)*mfexp(log_avg_F_GR+log_F_dev_end_GR);   
    
  F_end=F_end_L;
  F_end_apex=max(F_end);
  
  sel_wgted_tot=F_end/F_end_apex;
  sel_wgted_L=elem_prod(sel_wgted_tot, elem_div(F_end_L,F_end));
  
  wgt_wgted_L_denom=F_comm_prop+F_GR_prop;  //+F_HB_prop+F_cL_prop
  wgt_wgted_L_klb=F_comm_prop/wgt_wgted_L_denom*wholewgt_comm_klb(endyr)+ 
				  //F_cL_prop/wgt_wgted_L_denom*wholewgt_cL_klb(endyr)+ 
				  F_GR_prop/wgt_wgted_L_denom*wholewgt_GR_klb(endyr);                                      
  
FUNCTION get_msy
  
  //compute values as functions of F
  for(ff=1; ff<=n_iter_msy; ff++)
  {
    //uses fishery-weighted F's
    Z_age_msy=0.0;
    F_L_age_msy=0.0;
          
    F_L_age_msy=F_msy(ff)*sel_wgted_L; 
    Z_age_msy=M+F_L_age_msy;        
    
    N_age_msy(1)=1.0;
    for (iage=2; iage<=nages; iage++)
      {N_age_msy(iage)=N_age_msy(iage-1)*mfexp(-1.*Z_age_msy(iage-1));}
    N_age_msy(nages)=N_age_msy(nages)/(1.0-mfexp(-1.*Z_age_msy(nages)));
    N_age_msy_spawn(1,(nages-1))=elem_prod(N_age_msy(1,(nages-1)),
                                   mfexp((-1.*Z_age_msy(1,(nages-1)))*spawn_time_frac));                 
    N_age_msy_spawn(nages)=(N_age_msy_spawn(nages-1)*(mfexp(-1.*(Z_age_msy(nages-1)*(1.0-spawn_time_frac) + 
                            Z_age_msy(nages)*spawn_time_frac) )))/(1.0-mfexp(-1.*Z_age_msy(nages)));
                     
    spr_msy(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
	        
    //R_eq(ff)=SR_eq_func(R0, steep, spr_msy(1), spr_msy(ff), BiasCor, SR_switch);
    R_eq(ff)=BiasCor*R0;
	
    if (R_eq(ff)<dzero) {R_eq(ff)=dzero;}    
    N_age_msy*=R_eq(ff);
    N_age_msy_spawn*=R_eq(ff);
    
    for (iage=1; iage<=nages; iage++)
    {
      L_age_msy(iage)=N_age_msy(iage)*(F_L_age_msy(iage)/Z_age_msy(iage))*
                      (1.-mfexp(-1.*Z_age_msy(iage)));             
    }
    
    SSB_eq(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
	SSB_eq_knum(ff)=sum(elem_prod(N_age_msy_spawn,reprodknum));
	B_eq(ff)=sum(elem_prod(N_age_msy,wgt_mt));
    L_eq_klb(ff)=sum(elem_prod(L_age_msy,wgt_wgted_L_klb)); //in whole weight
    L_eq_knum(ff)=sum(L_age_msy)/1000.0;   
  }  
  
  msy_klb_out=max(L_eq_klb); //msy in whole weight 
  
  for(ff=1; ff<=n_iter_msy; ff++)
  {
   if(L_eq_klb(ff) == msy_klb_out) 
      {    
        SSB_msy_out=SSB_eq(ff);
        B_msy_out=B_eq(ff);
        R_msy_out=R_eq(ff);
        msy_knum_out=L_eq_knum(ff);      
        F_msy_out=F_msy(ff);  
        spr_msy_out=spr_msy(ff);      
      }
  }

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_per_recruit_stuff

  //static per-recruit stuff
 
  for(iyear=styr; iyear<=endyr; iyear++)
  {
    N_age_spr(1)=1.0;
    for(iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z(iyear,iage-1));}
    N_age_spr(nages)=N_age_spr(nages)/(1.0-mfexp(-1.*Z(iyear,nages)));    
    N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                mfexp(-1.*Z(iyear)(1,(nages-1))*spawn_time_frac));
    N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z(iyear)(nages-1)*(1.0-spawn_time_frac) + Z(iyear)(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z(iyear)(nages)));           
    spr_static(iyear)=sum(elem_prod(N_age_spr_spawn,reprod))/spr_F0;
  }
 
  //compute SSB/R and YPR as functions of F
  for(ff=1; ff<=n_iter_spr; ff++)
  {
    //uses fishery-weighted F's, same as in MSY calculations
    Z_age_spr=0.0;
    F_L_age_spr=0.0;

    F_L_age_spr=F_spr(ff)*sel_wgted_L;
    Z_age_spr=M+F_L_age_spr;

    N_age_spr(1)=1.0;
    for (iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));}
    N_age_spr(nages)=N_age_spr(nages)/(1-mfexp(-1.*Z_age_spr(nages)));
    N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));                 
    N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));
    spr_spr(ff)=sum(elem_prod(N_age_spr_spawn,reprod));
    L_spr(ff)=0.0;
    for (iage=1; iage<=nages; iage++)
    {
      L_age_spr(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      L_spr(ff)+=L_age_spr(iage)*wgt_wgted_L_klb(iage)*1000.0; //in lb whole wgt
    }   
  }
  spr_ratio=spr_spr/spr_F0;
  F20_dum=min(fabs(spr_ratio-0.2));
  F30_dum=min(fabs(spr_ratio-0.3));
  F40_dum=min(fabs(spr_ratio-0.4));
  for(ff=1; ff<=n_iter_spr; ff++)
  {   
      if (fabs(spr_ratio(ff)-0.2)==F20_dum) {F20_out=F_spr(ff);}	  
	  if (fabs(spr_ratio(ff)-0.3)==F30_dum) {F30_out=F_spr(ff);}	  
	  if (fabs(spr_ratio(ff)-0.4)==F40_dum) {F40_out=F_spr(ff);}
  }
  rec=column(N,1);
  rec_mean=sum(rec(styr_rec_spr, endyr_rec_spr))/nyrs_rec_spr;
  R_F30_out=rec_mean;
  F_L_age_spr=F30_out*sel_wgted_L;
  Z_age_spr=M+F_L_age_spr;

  N_age_spr(1)=R_F30_out;
  for (iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));}
  N_age_spr(nages)=N_age_spr(nages)/(1-mfexp(-1.*Z_age_spr(nages)));
  N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));                 
  N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));

  for (iage=1; iage<=nages; iage++)
    {
      L_age_F30(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));                 
    }
  SSB_F30_out=sum(elem_prod(N_age_spr_spawn,reprod));
  SSB_F30_knum_out=sum(elem_prod(N_age_spr_spawn,reprodknum));
  B_F30_out=sum(elem_prod(N_age_spr,wgt_mt));
  L_F30_klb_out=sum(elem_prod(L_age_F30,wgt_wgted_L_klb)); //in whole weight
  L_F30_knum_out=sum(L_age_F30)/1000.0; 
 

 //F40 calcs
  rec=column(N,1);
  rec_mean=sum(rec(styr_rec_spr, endyr_rec_spr))/nyrs_rec_spr;
  R_F40_out=rec_mean;
  F_L_age_spr=F40_out*sel_wgted_L;
  Z_age_spr=M+F_L_age_spr;

  N_age_spr(1)=R_F40_out;
  for (iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));}
  N_age_spr(nages)=N_age_spr(nages)/(1-mfexp(-1.*Z_age_spr(nages)));
  N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));                 
  N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));

  for (iage=1; iage<=nages; iage++)
    {
      L_age_F40(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));                 
    }
  SSB_F40_out=sum(elem_prod(N_age_spr_spawn,reprod));
  SSB_F40_knum_out=sum(elem_prod(N_age_spr_spawn,reprodknum));
  B_F40_out=sum(elem_prod(N_age_spr,wgt_mt));
  L_F40_klb_out=sum(elem_prod(L_age_F40,wgt_wgted_L_klb)); //in whole weight
  L_F40_knum_out=sum(L_age_F40)/1000.0; 
 
	
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_miscellaneous_stuff

//switch here if var_rec_dev <=dzero 
  if(var_rec_dev>0.0)
   {sigma_rec_dev=sqrt(var_rec_dev);} //sample SD of predicted residuals (may not equal rec_sigma)  
   else{sigma_rec_dev=0.0;}

  len_cv=elem_div(len_sd,meanlen_TL);
  len_cv_L=elem_div(len_sd_L,meanlen_TL_L);
  len_cv_F=elem_div(len_sd_F,meanlen_TL_F);
  
  //compute total landings- and discards-at-age in 1000 fish and klb whole weight
  L_total_num.initialize();
  L_total_klb.initialize();
  L_total_knum_yr.initialize();
  L_total_klb_yr.initialize();  
  
  for(iyear=styr; iyear<=endyr; iyear++)
  {
        L_total_klb_yr(iyear)=pred_comm_L_klb(iyear)+pred_GR_L_klb(iyear);//+pred_HB_L_klb(iyear)+pred_cL_L_klb(iyear)
        L_total_knum_yr(iyear)=pred_comm_L_knum(iyear)+pred_GR_L_knum(iyear);//+pred_HB_L_knum(iyear)+pred_cL_L_knum(iyear)
                
        B(iyear)=elem_prod(N(iyear),wgt_mt);
        totN(iyear)=sum(N(iyear));
        totB(iyear)=sum(B(iyear));                           
  }
  
  L_total_num=L_comm_num+L_GR_num;//+L_HB_num+L_cL_num   //landings at age in number fish
  L_total_klb=L_comm_klb+L_GR_klb;//+L_HB_klb+L_cL_klb   //landings at age in klb whole weight
 
  //Time series of interest
  B(endyr+1)=elem_prod(N(endyr+1),wgt_mt);
  totN(endyr+1)=sum(N(endyr+1));
  totB(endyr+1)=sum(B(endyr+1));  
  SdS0=SSB/S0;
  
  Fend_mean_temp=1.0;
  for (iyear=1; iyear<=selpar_n_yrs_wgted; iyear++) {Fend_mean_temp*=Fapex(endyr-iyear+1);}
  Fend_mean=pow(Fend_mean_temp,(1.0/selpar_n_yrs_wgted));	  
  if(F_msy_out>0)
    {
      FdF_msy=Fapex/F_msy_out;
      FdF_msy_end=FdF_msy(endyr);
      FdF_msy_end_mean=Fend_mean/F_msy_out;
    }
  if(SSB_msy_out>0)
    {
      SdSSB_msy=SSB/SSB_msy_out;
      SdSSB_msy_end=SdSSB_msy(endyr);
    }  

	if(F30_out>0)
    {
	  FdF30=Fapex/F30_out;
	  FdF30_end_mean=Fend_mean/F30_out;
	}
  if(SSB_F30_out>0)
    {
      SdSSB_F30=SSB/SSB_F30_out;
	  Sdmsst_F30=SSB/(smsy2msst75*SSB_F30_out);
      SdSSB_F30_end=SdSSB_F30(endyr);
	  Sdmsst_F30_end=Sdmsst_F30(endyr);
    }  	
	
	if(F40_out>0)
    {
	  FdF40=Fapex/F40_out;
	  FdF40_end_mean=Fend_mean/F40_out;
	}
  if(SSB_F40_out>0)
    {
      SdSSB_F40=SSB/SSB_F40_out;
	  Sdmsst_F40=SSB/(smsy2msst75*SSB_F40_out);
      SdSSB_F40_end=SdSSB_F40(endyr);
	  Sdmsst_F40_end=Sdmsst_F40(endyr);
    }  	
   //fill in log recruitment deviations for yrs they are nonzero
   for(iyear=styr_rec_dev; iyear<=endyr_rec_dev; iyear++)
     {log_rec_dev_output(iyear)=log_rec_dev(iyear);}
   //fill in log Nage deviations for ages they are nonzero (ages2+)
   for(iage=2; iage<=nages; iage++)
     {log_Nage_dev_output(iage)=log_Nage_dev(iage);}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_projection
  
    switch(Fproj_switch){
       case 1: //F=Fcurrent
          F_reg_proj=Fend_mean;
          break;
       case 2: //F=Fmsy
          F_reg_proj=F_msy_out;
          break;
       case 3: //F=F30
          F_reg_proj=F30_out;
          break;     
       case 4: //F=F40
          F_reg_proj=F40_out;
          break;          		  
       default: // no such switch available
          cout << "Error in input: Projection switch Fproj_switch must be set to 1, 2, 3, or 4." << endl;
          cout << "Presently it is set to " << Fproj_switch <<"."<< endl;
          exit(0);          
   }

  N_proj(styr_proj)=N(endyr+1); //initial conditions computed previously
 
  for (iyear=styr_proj; iyear<=endyr_proj; iyear++) //recruitment follows S-R curve (with bias correction) exactly
  {     
        if (iyear<styr_regs) {F_proj(iyear)=Fend_mean;}
		else {F_proj(iyear)=Fproj_mult*F_reg_proj;}
		
		FL_age_proj=sel_wgted_L*F_proj(iyear);
		
        Z_proj(iyear)=M+FL_age_proj;//+FD_age_proj;
        N_spawn_proj(iyear)(1,nages)=elem_prod(N_proj(iyear)(1,nages),(mfexp(-1.*(Z_proj(iyear)(1,nages))*spawn_time_frac))); //peak spawning time
		SSB_proj(iyear)= sum(elem_prod(N_spawn_proj(iyear),reprod));
        B_proj(iyear)=sum(elem_prod(N_proj(iyear),wgt_mt)); //uses spawning weight
	     
		for (iage=1; iage<=nages; iage++)
			{L_age_proj(iyear,iage)=N_proj(iyear,iage)*FL_age_proj(iage)*(1.-mfexp(-1.*Z_proj(iyear,iage)))/Z_proj(iyear,iage);
			}          
        L_knum_proj(iyear)=sum(L_age_proj(iyear))/1000.0;
	    L_klb_proj(iyear)=sum(elem_prod(L_age_proj(iyear),wgt_wgted_L_klb));     //in 1000 lb
		
		if (iyear<endyr_proj) {
			N_proj(iyear+1,1)=BiasCor*R0; //Changed to move away from an SR relationship
			//N_proj(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB_proj(iyear),SR_switch);
			N_proj(iyear+1)(2,nages)=++elem_prod(N_proj(iyear)(1,nages-1),(mfexp(-1.*Z_proj(iyear)(1,nages-1))));
			N_proj(iyear+1,nages)+=N_proj(iyear,nages)*mfexp(-1.*Z_proj(iyear,nages)); //plus group		
		}
  }
   R_proj=column(N_proj,1);                          
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   


FUNCTION evaluate_objective_function
  //fval=square(xdum-9.0);
  
  fval=0.0;
  fval_data=0.0;  
//---likelihoods---------------------------

//---Indices-------------------------------

  f_HB_cpue=0.0;
  f_HB_cpue=lk_lognormal(pred_HB_cpue, obs_HB_cpue, HB_cpue_cv, w_I_HB);
  fval+=f_HB_cpue;
  fval_data+=f_HB_cpue;   

//---Landings-------------------------------
  
  //f_comm_L in 1000 lb whole wgt  
  f_comm_L=lk_lognormal(pred_comm_L_klb(styr_comm_L,endyr_comm_L), obs_comm_L(styr_comm_L,endyr_comm_L),
                      comm_L_cv(styr_comm_L,endyr_comm_L), w_L);
  fval+=f_comm_L;
  fval_data+=f_comm_L;
  
  //f_GR_L in 1000 fish
  f_GR_L=lk_lognormal(pred_GR_L_knum(styr_GR_L,endyr_GR_L), obs_GR_L(styr_GR_L,endyr_GR_L), 
                      GR_L_cv(styr_GR_L,endyr_GR_L), w_L);
  fval+=f_GR_L;
  fval_data+=f_GR_L;  

//---Length comps-------------------------------

  //f_comm_lenc
  f_comm_lenc=lk_dirichlet_multinomial(nsamp_comm_lenc, pred_comm_lenc, obs_comm_lenc, nyr_comm_lenc, double(nlenbins), log_dm_comm_lc);
  fval+=f_comm_lenc;
  fval_data+=f_comm_lenc;
    
  
//---Age comps-------------------------------

//f_GR_agec
  f_GR_agec=lk_dirichlet_multinomial(nsamp_GR_agec, pred_GR_agec, obs_GR_agec, nyr_GR_agec, double(nages_agec), log_dm_GR_ac);
  fval+=f_GR_agec;
  fval_data+=f_GR_agec;

//-----------Constraints and penalties--------------------------------
  
  //Light penalty applied to log_Nage_dev for deviation from zero. If not estimated, this penalty equals zero.
  f_Nage_init=norm2(log_Nage_dev);        
  fval+=w_Nage_init*f_Nage_init;
  
  f_rec_dev=0.0;
  //rec_sigma_sq=square(rec_sigma);
  rec_logL_add=nyrs_rec*log(rec_sigma);
  f_rec_dev=(square(log_rec_dev(styr_rec_dev) + rec_sigma_sq/2.0)/(2.0*rec_sigma_sq));
  for(iyear=(styr_rec_dev+1); iyear<=endyr_rec_dev; iyear++)
  {f_rec_dev+=(square(log_rec_dev(iyear)-R_autocorr*log_rec_dev(iyear-1) + rec_sigma_sq/2.0)/
               (2.0*rec_sigma_sq));}
  f_rec_dev+=rec_logL_add;            
  fval+=w_rec*f_rec_dev;
   
  f_rec_dev_early=0.0; //possible extra constraint on early rec deviations
  if (w_rec_early>0.0)
    { if (styr_rec_dev<endyr_rec_phase1)
        {  
          for(iyear=styr_rec_dev; iyear<=endyr_rec_phase1; iyear++)
          //{f_rec_dev_early+=(square(log_rec_dev(iyear)-R_autocorr*log_rec_dev(iyear-1) + rec_sigma_sq/2.0)/
          //                  (2.0*rec_sigma_sq)) + rec_logL_add;}
          {f_rec_dev_early+=square(log_rec_dev(iyear));}
        }
  fval+=w_rec_early*f_rec_dev_early;
  }
  
  f_rec_dev_end=0.0; //possible extra constraint on ending rec deviations
  if (w_rec_end>0.0)
  { if (endyr_rec_phase2<endyr_rec_dev)
        {  
          for(iyear=(endyr_rec_phase2+1); iyear<=endyr_rec_dev; iyear++)
          //{f_rec_dev_end+=(square(log_rec_dev(iyear)-R_autocorr*log_rec_dev(iyear-1) + rec_sigma_sq/2.0)/
          //                 (2.0*rec_sigma_sq)) + rec_logL_add;}
          {f_rec_dev_end+=square(log_rec_dev(iyear));}
        }
      fval+=w_rec_end*f_rec_dev_end;
   }  

  //Ftune penalty: does not apply in last phase
  f_Ftune=0.0; 
  if (w_Ftune>0.0)
  {if (set_Ftune>0.0 && !last_phase()) {f_Ftune=square(Fapex(set_Ftune_yr)-set_Ftune);}
   fval+=w_Ftune*f_Ftune;
  }

  //Penalty if apical F exceeds 3.0
  f_fullF_constraint=0.0;
  if (w_fullF>0.0)
  {for (iyear=styr; iyear<=endyr; iyear++)
        {if(Fapex(iyear)>3.0) {f_fullF_constraint+=(mfexp(Fapex(iyear)-3.0)-1.0);}}
   fval+=w_fullF*f_fullF_constraint;
  }
  
 //Random walk components of fishery dependent indices
 //f_comm_RWq_cpue=0.0;
 //for (iyear=styr_comm_cpue; iyear<endyr_comm_cpue; iyear++)
 //    {f_comm_RWq_cpue+=square(q_RW_log_dev_comm(iyear))/(2.0*set_RWq_var);}
 //fval+=f_comm_RWq_cpue;
 //
 //f_cL_RWq_cpue=0.0;
 //for (iyear=styr_cL_cpue; iyear<endyr_cL_cpue; iyear++)
 //    {f_cL_RWq_cpue+=square(q_RW_log_dev_cL(iyear))/(2.0*set_RWq_var);}
 //fval+=f_cL_RWq_cpue; 
 //
 f_HB_RWq_cpue=0.0;
 for (iyear=styr_HB_cpue; iyear<endyr_HB_cpue; iyear++)
     {f_HB_RWq_cpue+=square(q_RW_log_dev_HB(iyear))/(2.0*set_RWq_var);}
 fval+=f_HB_RWq_cpue;     
  
//---Priors---------------------------------------------------
//neg_log_prior arguments: estimate, prior mean, prior var/-CV, pdf type
//Variance input as a negative value is considered to be CV in arithmetic space (CV=-1 implies loose prior) 
//pdf type 1=none, 2=lognormal, 3=normal, 4=beta 
  f_priors=0.0; 
  f_priors+=neg_log_prior(len_cv_val,set_len_cv(5),set_len_cv(6),set_len_cv(7));
   
  f_priors+=neg_log_prior(steep,set_steep(5),set_steep(6),set_steep(7)); 
  f_priors+=neg_log_prior(log_R0,set_log_R0(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(R_autocorr,set_R_autocorr(5),set_R_autocorr(6),set_R_autocorr(7));
  f_priors+=neg_log_prior(rec_sigma,set_rec_sigma(5),set_rec_sigma(6),set_rec_sigma(7));
 
  f_priors+=neg_log_prior(selpar_A50_comm1,set_selpar_A50_comm1(5), set_selpar_A50_comm1(6), set_selpar_A50_comm1(7));
  f_priors+=neg_log_prior(selpar_slope_comm1,set_selpar_slope_comm1(5), set_selpar_slope_comm1(6), set_selpar_slope_comm1(7));
  //f_priors+=neg_log_prior(selpar_A50_comm2,set_selpar_A50_comm2(5), set_selpar_A50_comm2(6), set_selpar_A50_comm2(7));
  //f_priors+=neg_log_prior(selpar_slope_comm2,set_selpar_slope_comm2(5), set_selpar_slope_comm2(6), set_selpar_slope_comm2(7));
  // f_priors+=neg_log_prior(selpar_A502_comm2,set_selpar_A502_comm2(5), set_selpar_A502_comm2(6), set_selpar_A502_comm2(7));
  // f_priors+=neg_log_prior(selpar_slope2_comm2,set_selpar_slope2_comm2(5), set_selpar_slope2_comm2(6), set_selpar_slope2_comm2(7));  
  //f_priors+=neg_log_prior(selpar_A50_comm3,set_selpar_A50_comm3(5), set_selpar_A50_comm3(6), set_selpar_A50_comm3(7));
  //f_priors+=neg_log_prior(selpar_slope_comm3,set_selpar_slope_comm3(5), set_selpar_slope_comm3(6), set_selpar_slope_comm3(7));  


  f_priors+=neg_log_prior(selpar_A50_GR1,set_selpar_A50_GR1(5), set_selpar_A50_GR1(6), set_selpar_A50_GR1(7));
  f_priors+=neg_log_prior(selpar_slope_GR1,set_selpar_slope_GR1(5), set_selpar_slope_GR1(6), set_selpar_slope_GR1(7));
  f_priors+=neg_log_prior(selpar_A50_GR2,set_selpar_A50_GR2(5), set_selpar_A50_GR2(6), set_selpar_A50_GR2(7));
  f_priors+=neg_log_prior(selpar_slope_GR2,set_selpar_slope_GR2(5), set_selpar_slope_GR2(6), set_selpar_slope_GR2(7));
  
   f_priors+=neg_log_prior(log_q_HB,set_log_q_HB(5),set_log_q_HB(6),set_log_q_HB(7));
      
  f_priors+=neg_log_prior(log_dm_comm_lc,set_log_dm_comm_lc(5),set_log_dm_comm_lc(6),set_log_dm_comm_lc(7));
  //f_priors+=neg_log_prior(log_dm_cL_lc,set_log_dm_cL_lc(5),set_log_dm_cL_lc(6),set_log_dm_cL_lc(7));
  f_priors+=neg_log_prior(log_dm_GR_ac,set_log_dm_GR_ac(5),set_log_dm_GR_ac(6),set_log_dm_GR_ac(7));
  //f_priors+=neg_log_prior(log_dm_GR_lc,set_log_dm_GR_lc(5),set_log_dm_GR_lc(6),set_log_dm_GR_lc(7));
  
  f_priors+=neg_log_prior(F_init,set_F_init(5),set_F_init(6),set_F_init(7));  
  
  fval+=f_priors;

  //Indv Liklihood Components (for WAIC calcs)
  for (int i=styr_HB_cpue;i<=endyr_HB_cpue;i++){
   LL_HB_cpue(i)=log(obs_HB_cpue(i)+0.0001)+0.5*log(2*3.141593)+log(sqrt(log(1.0+square(HB_cpue_cv(i)/w_I_HB))))+square(log((obs_HB_cpue(i)+0.0001)/(pred_HB_cpue(i)+0.0001)))/(2*log(1.0+square(HB_cpue_cv(i)/w_I_HB)));
  }
  for (int i=styr_comm_L;i<=endyr_comm_L;i++){
   LL_comm_L(i)=log(obs_comm_L(i)+0.0001)+0.5*log(2*3.141593)+log(sqrt(log(1.0+square(comm_L_cv(i)/w_L))))+square(log((obs_comm_L(i)+0.0001)/(pred_comm_L_klb(i)+0.0001)))/(2*log(1.0+square(comm_L_cv(i)/w_L)));
  }
  for (int i=styr_GR_L;i<=endyr_GR_L;i++){
   LL_GR_L(i)=log(obs_GR_L(i)+0.0001)+0.5*log(2*3.141593)+log(sqrt(log(1.0+square(GR_L_cv(i)/w_L))))+square(log((obs_GR_L(i)+0.0001)/(pred_GR_L_knum(i)+0.0001)))/(2*log(1.0+square(GR_L_cv(i)/w_L)));
  }

  LL_comm_lc=lk_dirichlet_multinomial(nsamp_comm_lenc, pred_comm_lenc, obs_comm_lenc, nyr_comm_lenc, double(nlenbins), log_dm_comm_lc);

  for (int i=1;i<=26;i++){
   LL_GR_ac(i)=lk_dirichlet_multinomial_ind(nsamp_GR_agec(i), pred_GR_agec(i), obs_GR_agec(i), double(nages_agec), log_dm_GR_ac);
  }

//----------------------------------------------------------------------------------
//Logistic function: 2 parameters
FUNCTION dvar_vector logistic(const dvar_vector& ages, const dvariable& A50, const dvariable& slope)
  //ages=vector of ages, A50=age at 50% selectivity, slope=rate of increase
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1./(1.+mfexp(-1.*slope*(ages-A50))); //logistic;  
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;

//-----------------------------------------------------------------------------------                                                                           
//Logistic-exponential: 4 parameters (but 1 is fixed)
FUNCTION dvar_vector logistic_exponential(const dvar_vector& ages, const dvariable& A50, const dvariable& slope, const dvariable& sigma, const dvariable& joint)
  //ages=vector of ages, A50=age at 50% sel (ascending limb), slope=rate of increase, sigma=controls rate of descent (descending)                               
  //joint=age to join curves                                                                                                                                    
  RETURN_ARRAYS_INCREMENT();                                                                                                                                    
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());                                                                                                         
  Sel_Tmp=1.0;                                                                                                                                                  
  for (iage=1; iage<=nages; iage++)                                                                                                                             
  {                                                                                                                                                             
   if (ages(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope*(ages(iage)-A50)));}                                                                             
   if (ages(iage)>joint){Sel_Tmp(iage)=mfexp(-1.*square((ages(iage)-joint)/sigma));}                                                                            
  }                                                                                                                                                             
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);                                                                                                                                 
  RETURN_ARRAYS_DECREMENT();                                                                                                                                    
  return Sel_Tmp;   

//-----------------------------------------------------------------------------------
//Logistic function: 4 parameters
FUNCTION dvar_vector logistic_double(const dvar_vector& ages, const dvariable& A501, const dvariable& slope1, const dvariable& A502, const dvariable& slope2)
  //ages=vector of ages, A50=age at 50% selectivity, slope=rate of increase, A502=age at 50% decrease additive to A501, slope2=slope of decrease
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=elem_prod( (1./(1.+mfexp(-1.*slope1*(ages-A501)))),(1.-(1./(1.+mfexp(-1.*slope2*(ages-(A501+A502)))))) );     
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;

//-----------------------------------------------------------------------------------
//Jointed logistic function: 6 parameters (increasing and decreasing logistics joined at peak selectivity)
FUNCTION dvar_vector logistic_joint(const dvar_vector& ages, const dvariable& A501, const dvariable& slope1, const dvariable& A502, const dvariable& slope2, const dvariable& satval, const dvariable& joint)
  //ages=vector of ages, A501=age at 50% sel (ascending limb), slope1=rate of increase,A502=age at 50% sel (descending), slope1=rate of increase (ascending), 
  //satval=saturation value of descending limb, joint=location in age vector to join curves (may equal age or age + 1 if age-0 is included)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1.0; 
  for (iage=1; iage<=nages; iage++)
  {
   if (double(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope1*(ages(iage)-A501)));}  
   if (double(iage)>joint){Sel_Tmp(iage)=1.0-(1.0-satval)/(1.+mfexp(-1.*slope2*(ages(iage)-A502)));}  
  }  
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;

//-----------------------------------------------------------------------------------  
//Double Gaussian function: 6 parameters (as in SS3)
FUNCTION dvar_vector gaussian_double(const dvar_vector& ages, const dvariable& peak, const dvariable& top, const dvariable& ascwid, const dvariable& deswid, const dvariable& init, const dvariable& final)
  //ages=vector of ages, peak=ascending inflection location (as logistic), top=width of plateau, ascwid=ascent width (as log(width))
  //deswid=descent width (as log(width))
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step1(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step2(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step3(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step4(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step5(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step6(ages.indexmin(),ages.indexmax());
  dvar_vector pars_tmp(1,6); dvar_vector sel_tmp_iq(1,2);
  
  pars_tmp(1)=peak;
  pars_tmp(2)=peak+1.0+(0.99*ages(nages)-peak-1.0)/(1.0+mfexp(-top));
  pars_tmp(3)=mfexp(ascwid);
  pars_tmp(4)=mfexp(deswid);
  pars_tmp(5)=1.0/(1.0+mfexp(-init));
  pars_tmp(6)=1.0/(1.0+mfexp(-final));
       
  sel_tmp_iq(1)=mfexp(-(square(ages(1)-pars_tmp(1))/pars_tmp(3)));
  sel_tmp_iq(2)=mfexp(-(square(ages(nages)-pars_tmp(2))/pars_tmp(4)));
  
  sel_step1=mfexp(-(square(ages-pars_tmp(1))/pars_tmp(3)));
  sel_step2=pars_tmp(5)+(1.0-pars_tmp(5))*(sel_step1-sel_tmp_iq(1))/(1.0-sel_tmp_iq(1));  
  sel_step3=mfexp(-(square(ages-pars_tmp(2))/pars_tmp(4)));
  sel_step4=1.0+(pars_tmp(6)-1.0)*(sel_step3-1.0)/(sel_tmp_iq(2)-1.0);
  sel_step5=1.0/ (1.0+mfexp(-(20.0* elem_div((ages-pars_tmp(1)), (1.0+sfabs(ages-pars_tmp(1)))) )));
  sel_step6=1.0/(1.0+mfexp(-(20.0*elem_div((ages-pars_tmp(2)),(1.0+sfabs(ages-pars_tmp(2)))) )));  

  Sel_Tmp=elem_prod(sel_step2,(1.0-sel_step5))+ 
          elem_prod(sel_step5,((1.0-sel_step6)+ elem_prod(sel_step4,sel_step6)) ); 
 
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
    
//-----------------------------------------------------------------------------------    
//Spawner-recruit function (Beverton-Holt or Ricker)
FUNCTION dvariable SR_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& SSB, int func)
  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, SSB=spawning biomass
  //func=1 for Beverton-Holt, 2 for Ricker
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  switch(func) {
    case 1: //Beverton-Holt
      Recruits_Tmp=((0.8*R0*h*SSB)/(0.2*R0*spr_F0*(1.0-h)+(h-0.2)*SSB));       
    break;
    case 2: //Ricker
      Recruits_Tmp=((SSB/spr_F0)*mfexp(h*(1-SSB/(R0*spr_F0))));       
    break;
	case 3: //None
      Recruits_Tmp=R0;      
    break;
  }
  RETURN_ARRAYS_DECREMENT();
  return Recruits_Tmp;
  
//-----------------------------------------------------------------------------------    
//Spawner-recruit equilibrium function (Beverton-Holt or Ricker)
FUNCTION dvariable SR_eq_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& spr_F, const dvariable& BC, int func)
  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, spr_F=spawners per recruit @ F, BC=bias correction
  //func=1 for Beverton-Holt, 2 for Ricker
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  switch(func) {
    case 1: //Beverton-Holt
      Recruits_Tmp=(R0/((5.0*h-1.0)*spr_F))*(BC*4.0*h*spr_F-spr_F0*(1.0-h));    
    break;
    case 2: //Ricker
      Recruits_Tmp=R0/(spr_F/spr_F0)*(1.0+log(BC*spr_F/spr_F0)/h);      
    break;
	case 3: //None
      Recruits_Tmp=BC*R0;      
    break;
  }
  RETURN_ARRAYS_DECREMENT();
  return Recruits_Tmp;

//-----------------------------------------------------------------------------------
//compute multinomial effective sample size for a single yr
FUNCTION dvariable multinom_eff_N(const dvar_vector& pred_comp, const dvar_vector& obs_comp)
  //pred_comp=vector of predicted comps, obscomp=vector of observed comps
  dvariable EffN_Tmp; dvariable numer; dvariable denom;
  RETURN_ARRAYS_INCREMENT();
  numer=sum( elem_prod(pred_comp,(1.0-pred_comp)) );
  denom=sum( square(obs_comp-pred_comp) );
  if (denom>0.0) {EffN_Tmp=numer/denom;}
  else {EffN_Tmp=-missing;}                            
  RETURN_ARRAYS_DECREMENT();
  return EffN_Tmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: lognormal
FUNCTION dvariable lk_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
  //pred=vector of predicted vals, obs=vector of observed vals, cv=vector of CVs in arithmetic space, wgt_dat=constant scaling of CVs
  //small_number is small value to avoid log(0) during search
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.0001;
  dvariable pi=3.141593;
  dvariable n=(obs.indexmax()-obs.indexmin())+1;
  dvar_vector var(cv.indexmin(),cv.indexmax()); //variance in log space
  var=log(1.0+square(cv/wgt_dat));   // convert cv in arithmetic space to variance in log space
  //LkvalTmp=sum(0.5*elem_div(square(log(elem_div((pred+small_number),(obs+small_number)))),var) );

  ///*
  //Full Likelihood
  LkvalTmp=sum(log(obs+small_number))+n*0.5*log(2*pi)+sum(log(sqrt(var)))+sum(elem_div(square(log(elem_div((obs+small_number),(pred+small_number)))),2*var));
  /*
  LkvalTmp=0.0;
  for (int ii=obs.indexmin(); ii<=obs.indexmax(); ii++){
   LkvalTmp += log(obs(ii)+small_number)+0.5*log(2*pi)+log(sqrt(var(ii)))+square(log(obs(ii)+small_number)-log(pred(ii)+small_number))/(2*var(ii));
  }
  */
  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: multinomial
FUNCTION dvariable lk_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const double& minSS, const dvariable& wgt_dat)
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.0001;
  LkvalTmp=0.0;
  /*
  for (int ii=1; ii<=ncomp; ii++){
   if (nsamp(ii)>=minSS){
    LkvalTmp-=wgt_dat*nsamp(ii)*sum(elem_prod((obs_comp(ii)+small_number),log(elem_div((pred_comp(ii)+small_number),(obs_comp(ii)+small_number)))));
    }
  }
  */

  //Multinomial Likelihood
  for (int ii=1; ii<=ncomp; ii++){
   if (nsamp(ii)>=minSS){
   LkvalTmp -= wgt_dat*nsamp(ii)*sum(elem_prod((obs_comp(ii)+small_number),log(pred_comp(ii)+small_number)));
   }
  }
  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: robust multinomial
FUNCTION dvariable lk_robust_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& wgt_dat)
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.0001;
  LkvalTmp=0.0;
  dvar_matrix Eprime=elem_prod((1.0-obs_comp), obs_comp)+0.1/mbin; //E' of Francis 2011, p.1131  
  dvar_vector nsamp_wgt=nsamp*wgt_dat;
  //cout<<nsamp_wgt<<endl;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {LkvalTmp+= sum(0.5*log(Eprime(ii))-log(small_number+mfexp(elem_div((-square(obs_comp(ii)-pred_comp(ii))) , (Eprime(ii)*2.0/nsamp_wgt(ii)) ))) );
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;
  
//-----------------------------------------------------------------------------------
//Likelihood contribution: Dirichlet-multinomial, Linear Parameterization
FUNCTION dvariable lk_dirichlet_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const double& mbin, const dvariable& log_dir_par)
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.0001;
  LkvalTmp=0.0;
  dvar_matrix obs_comp_plus=obs_comp;
  dvar_matrix obs_comp_norm=obs_comp*0.0;
  dvar_vector Beta=nsamp*mfexp(log_dir_par);   //number of samples * the DM param
  dvar_vector ones(1,mbin);

  for (int jj=1;jj<=mbin;jj++){
   ones(jj) = 1;   //Filling in ones vector
  }

  for (int ii=1;ii<=ncomp;ii++){
    for (int jj=1;jj<=mbin;jj++){
     if(obs_comp(ii,jj)==0){
      obs_comp_plus(ii,jj) = small_number;  //Suppressing zeroes in the observed data. Adding very small number to zero observations, then renormalizing
     }
    }
   obs_comp_norm(ii) = obs_comp_plus(ii)/sum(obs_comp_plus(ii));
  }
  
   //dvar_vector nsamp_adjust=mfexp(log_dir_par);
  for (int ii=1; ii<=ncomp; ii++)
  {
   //Full Likelihood
   LkvalTmp += -1*(gammln(nsamp(ii)+1)-sum(gammln(nsamp(ii)*obs_comp_norm(ii)+ones))+gammln(Beta(ii))-gammln(nsamp(ii)+Beta(ii))+sum(gammln(nsamp(ii)*obs_comp_norm(ii)+Beta(ii)*pred_comp(ii))-gammln(Beta(ii)*pred_comp(ii))));
  }
  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

//-----------------------------------------------------------------------------------
//Individual Likelihood contribution: Dirichlet-multinomial, Linear Parameterization
FUNCTION dvariable lk_dirichlet_multinomial_ind(const double& nsamp, const dvar_vector& pred_comp, const dvar_vector& obs_comp, const double& mbin, const dvariable& log_dir_par)
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.0001;
  LkvalTmp=0.0;
  dvar_vector obs_comp_plus=obs_comp;
  dvar_vector obs_comp_norm=obs_comp*0.0;
  dvariable Beta=nsamp*mfexp(log_dir_par);   //number of samples * the DM param
  dvar_vector ones(1,mbin);

  for (int jj=1;jj<=mbin;jj++){
   ones(jj) = 1;   //Filling in ones vector
  }

    for (int jj=1;jj<=mbin;jj++){
     if(obs_comp(jj)==0){
      obs_comp_plus(jj) = small_number;  //Suppressing zeroes in the observed data. Adding very small number to zero observations, then renormalizing
     }
    }
    
   obs_comp_norm = obs_comp_plus/sum(obs_comp_plus);
  
   //Full Likelihood
   LkvalTmp = -1.0*(gammln(nsamp+1)-sum(gammln(nsamp*obs_comp_norm+ones))+gammln(Beta)-gammln(nsamp+Beta)+sum(gammln(nsamp*obs_comp_norm+Beta*pred_comp)-gammln(Beta*pred_comp)));
  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: logistic normal (aka multivariate logistic in iSCAM; logistic normal in Francis' terminology)
FUNCTION dvariable lk_logistic_normal(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS)
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.0001;
  LkvalTmp=0.0;
  dvar_matrix nu=pred_comp+0.0;
  dvar_matrix pred_plus=pred_comp+small_number;
  dvar_matrix obs_plus=obs_comp+small_number;

  dvariable nu_mean;
  dvariable nu_sum_sq;
  dvariable tau_hat_sq;
  dvariable year_count; //keeps track of years included in likelihood (i.e., that meet the sample size requirement)
   
  LkvalTmp=0.0;
  nu_sum_sq=0.0;
  year_count=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {
		year_count+=1.0;
		nu_mean=sum( log(obs_plus(ii))-log(pred_plus(ii))  )/mbin;	//year-specific mean log residual
		for (int jj=1; jj<=mbin;jj++)
		{
			nu(ii,jj) = log(obs_plus(ii,jj)) - log(pred_plus(ii,jj)) - nu_mean;
			nu_sum_sq += square(nu(ii,jj));
		}
    }
  }  
  if (year_count>0.0)
  {
	  tau_hat_sq = nu_sum_sq/((mbin-1.0)*year_count);
	  LkvalTmp = (mbin-1.0)*year_count*log(tau_hat_sq);
  }
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;
  
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//Likelihood contribution: priors
FUNCTION  dvariable neg_log_prior(dvariable pred, const double& prior, dvariable var, int pdf)
  //prior=prior point estimate, var=variance (if negative, treated as CV in arithmetic space), pred=predicted value, pdf=prior type (1=none, 2=lognormal, 3=normal, 4=beta)
    dvariable LkvalTmp;
    dvariable alpha, beta, ab_iq;
    dvariable big_number=1e10;
    dvariable pi=3.141593;
    LkvalTmp=0.0;
    // compute generic pdf's
    switch(pdf) {
        case 1: //option to turn off prior
          LkvalTmp=0.0;
          break;
        case 2: // lognormal
          if(prior<=0.0) cout << "YIKES: Don't use a lognormal distn for a negative prior" << endl;
          else if(pred<=0) LkvalTmp=big_number=1e10;
          else {
            if(var<0.0) var=log(1.0+var*var) ;      // convert cv to variance on log scale
//            LkvalTmp= 0.5*( square(log(pred/prior))/var + log(var) );
            LkvalTmp= 0.5*(square(log(prior/pred))/var) + log(sqrt(var)) + log(prior) + 0.5*log(2*pi) ; //Full Lognormal likelihood (including constants)
          }
        break;
        case 3: // normal
          if(var<0.0 && prior!=0.0) var=square(var*prior);       // convert cv to variance on observation scale
          else if(var<0.0 && prior==0.0) var=-var;               // cv not really appropriate if prior value equals zero
//          LkvalTmp= 0.5*( square(pred-prior)/var + log(var) );
          LkvalTmp= 0.5*(square(prior-pred)/var) + log(sqrt(var)) + 0.5*log(2*pi) ;   //Full Normal Likelihood
          break;
        case 4: // beta
          if(var<0.0) var=square(var*prior);          // convert cv to variance on observation scale
          if(prior<=0.0 || prior>=1.0) cout << "YIKES: Don't use a beta distn for a prior outside (0,1)" << endl;
          ab_iq=prior*(1.0-prior)/var - 1.0; alpha=prior*ab_iq; beta=(1.0-prior)*ab_iq;
          if(pred>=0 && pred<=1) LkvalTmp= (1.0-alpha)*log(pred)+(1.0-beta)*log(1.0-pred)-gammln(alpha+beta)+gammln(alpha)+gammln(beta);
          else LkvalTmp=big_number;
          break;
        default: // no such prior pdf currently available
          cout << "The prior must be either 1(lognormal), 2(normal), or 3(beta)." << endl;
          cout << "Presently it is " << pdf << endl;
          exit(0);
    }
    return LkvalTmp;

//-----------------------------------------------------------------------------------
//SDNR: age comp likelihood (assumes fits are done with the robust multinomial function)
FUNCTION dvariable sdnr_multinomial(const double& ncomp, const dvar_vector& ages, const dvar_vector& nsamp, 
                                    const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const dvariable& wgt_dat)
  //ncomp=number of years of data, ages=vector of ages, nsamp=vector of N's, 
  //pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, wgt_dat=likelihood weight for data source
  RETURN_ARRAYS_INCREMENT();
  dvariable SdnrTmp;
  dvar_vector o(1,ncomp);  
  dvar_vector p(1,ncomp);  
  dvar_vector ose(1,ncomp);  
  dvar_vector res(1,ncomp);
  SdnrTmp=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {
    o(ii)=sum(elem_prod(ages,obs_comp(ii)));
    p(ii)=sum(elem_prod(ages,pred_comp(ii)));
    ose(ii)=sqrt((sum(elem_prod(square(ages),pred_comp(ii)))-square(p(ii)))/(nsamp(ii)*wgt_dat));
  }
  res=elem_div((o-p),ose); 
  SdnrTmp=sqrt(sum(square(res-(sum(res)/ncomp))/(ncomp-1.0))); 
  RETURN_ARRAYS_DECREMENT();
  return SdnrTmp;

//-----------------------------------------------------------------------------------
//SDNR: lognormal likelihood
FUNCTION dvariable sdnr_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
  //nyr=number of years of data, pred=vector of predicted data, obs=vector of observed data, cv=vector of cv's, wgt_dat=likelihood weight for data source
  RETURN_ARRAYS_INCREMENT();
  dvariable SdnrTmp;
  dvariable small_number=0.00001;
  dvariable n;
  dvar_vector res(cv.indexmin(),cv.indexmax());
  SdnrTmp=0.0;
  res=elem_div(log(elem_div(obs+small_number,pred+small_number)),sqrt(log(1+square(cv/wgt_dat))));
  n=cv.indexmax()-cv.indexmin()+1;
  SdnrTmp=sqrt(sum(square(res-(sum(res)/n))/(n-1.0))); 
  RETURN_ARRAYS_DECREMENT();
  return SdnrTmp; 

//-----------------------------------------------------------------------------------
REPORT_SECTION

  if (last_phase())  
  {
      cout<<"start report"<<endl;
       //cout<<"xdum = "<<xdum<<endl;
	  get_weighted_current();
       cout<<"got weighted"<<endl;
      get_msy();
       cout<<"got msy"<<endl;
      get_per_recruit_stuff();
       cout<<"got per recruit"<<endl;  
      get_miscellaneous_stuff();
       cout<<"got misc stuff"<<endl;
	  get_projection();
	   cout<<"got projection"<<endl;
	  
      grad_max=objective_function_value::pobjfun->gmax;
      time(&finish);
	  elapsed_time=difftime(finish,start);
	  hour=long(elapsed_time)/3600;
	  minute=long(elapsed_time)%3600/60;
	  second=(long(elapsed_time)%3600)%60;
	  cout<<endl<<endl<<"*******************************************"<<endl;
	  cout<<"--Start time: "<<ctime(&start)<<endl;
	  cout<<"--Finish time: "<<ctime(&finish)<<endl;
	  cout<<"--Runtime: ";
	  cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	  cout << "--TotalLikelihood: " << fval << endl;
      cout<<"--Final gradient: "<<objective_function_value::pobjfun->gmax << endl;
	  cout<<"*******************************************"<<endl;
      
      cout <<endl;     
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;
      cout << "F40=" << F40_out<< "   SSB.F40=" << SSB_F40_out <<endl;
      cout <<"F status="<<FdF_msy_end<<endl;
      cout <<"Pop status="<<SdSSB_msy_end<<endl;
      cout << "R0="<<R0<<endl;
      //cout << "xdum " << xdum << endl;
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;  
      
      report << "TotalLikelihood " << fval << endl;
      report << "N" << endl;
      report << N<<endl;
      report << "F" << endl;
      report << F <<endl;
      // report << "prob_belowsizelim_block3" << endl;
	  // report<<prob_belowsizelim_block3<<endl;
	        
      sdnr_I_HB=sdnr_lognormal(pred_HB_cpue, obs_HB_cpue, HB_cpue_cv, w_I_HB);      
      
      //#################################################################################################
      //##  Passing parameters to vector for bounds check plotting
      //################################################################################################# 
       Linf_out(8)=Linf; Linf_out(1,7)=set_Linf; 
       K_out(8)=K; K_out(1,7)=set_K;
       t0_out(8)=t0; t0_out(1,7)=set_t0;
       len_cv_val_out(8)=len_cv_val; len_cv_val_out(1,7)=set_len_cv;
	   
	   Linf_L_out(8)=Linf_L; Linf_L_out(1,7)=set_Linf_L; 
       K_L_out(8)=K_L; K_L_out(1,7)=set_K_L;
       t0_L_out(8)=t0_L; t0_L_out(1,7)=set_t0_L;
       len_cv_val_L_out(8)=len_cv_val_L; len_cv_val_L_out(1,7)=set_len_cv_L;
	   
	   Linf_F_out(8)=Linf_F; Linf_F_out(1,7)=set_Linf_F; 
       K_F_out(8)=K_F; K_F_out(1,7)=set_K_F;
       t0_F_out(8)=t0_F; t0_F_out(1,7)=set_t0_F;
       len_cv_val_F_out(8)=len_cv_val_F; len_cv_val_F_out(1,7)=set_len_cv_F;
	   
       log_R0_out(8)=log_R0; log_R0_out(1,7)=set_log_R0;
       M_constant_out(8)=M_constant; M_constant_out(1,7)=set_M_constant;
       steep_out(8)=steep; steep_out(1,7)=set_steep;
       rec_sigma_out(8)=rec_sigma; rec_sigma_out(1,7)=set_rec_sigma;
       R_autocorr_out(8)=R_autocorr; R_autocorr_out(1,7)=set_R_autocorr;
	   
	   log_dm_comm_lc_out(8)=log_dm_comm_lc; log_dm_comm_lc_out(1,7)=set_log_dm_comm_lc;
	   //log_dm_cL_lc_out(8)=log_dm_cL_lc; log_dm_cL_lc_out(1,7)=set_log_dm_cL_lc;
	   //log_dm_GR_lc_out(8)=log_dm_GR_lc; log_dm_GR_lc_out(1,7)=set_log_dm_GR_lc;
       log_dm_GR_ac_out(8)=log_dm_GR_ac; log_dm_GR_ac_out(1,7)=set_log_dm_GR_ac;
	   
       selpar_A50_comm1_out(8)=selpar_A50_comm1; selpar_A50_comm1_out(1,7)=set_selpar_A50_comm1;
       selpar_slope_comm1_out(8)=selpar_slope_comm1; selpar_slope_comm1_out(1,7)=set_selpar_slope_comm1;
       
	   selpar_A50_GR1_out(8)=selpar_A50_GR1; selpar_A50_GR1_out(1,7)=set_selpar_A50_GR1;
       selpar_slope_GR1_out(8)=selpar_slope_GR1; selpar_slope_GR1_out(1,7)=set_selpar_slope_GR1;
	   selpar_A50_GR2_out(8)=selpar_A50_GR2; selpar_A50_GR2_out(1,7)=set_selpar_A50_GR2;
       selpar_slope_GR2_out(8)=selpar_slope_GR2; selpar_slope_GR2_out(1,7)=set_selpar_slope_GR2;
        
       log_q_HB_out(8)=log_q_HB; log_q_HB_out(1,7)=set_log_q_HB;
                     
       log_avg_F_comm_out(8)=log_avg_F_comm; log_avg_F_comm_out(1,7)=set_log_avg_F_comm;
	   //log_avg_F_cL_out(8)=log_avg_F_cL; log_avg_F_cL_out(1,7)=set_log_avg_F_cL;
       log_avg_F_GR_out(8)=log_avg_F_GR; log_avg_F_GR_out(1,7)=set_log_avg_F_GR;
       F_init_out(8)=F_init; F_init_out(1,7)=set_F_init;       
        
       log_rec_dev_out(styr_rec_dev, endyr_rec_dev)=log_rec_dev;
       log_F_dev_comm_out(styr_comm_L,endyr_comm_L)=log_F_dev_comm;
	   //log_F_dev_cL_out(styr_cL_L,endyr_cL_L)=log_F_dev_cL;
       log_F_dev_GR_out(styr_GR_L,endyr_GR_L)=log_F_dev_GR;
           
   #include "co22_make_Robject4.cxx"   // write the R-compatible report

  } //endl last phase loop     
  
 
