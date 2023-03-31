//Recreation of Pacific Hake Assessment Model by N. Fisch July 2021

  //With Dirichlet-Multinomial Saturating for composition data 

DATA_SECTION  
  init_int age_err                              //ageing error indicator (0=off, 1=on)
  init_int b_ramp                               //bias ramp indicator (0=off, 1=on)

  init_int fyear                                //first year
  init_int lyear                                //last year
  init_int fage                                 //first age
  init_int lage                                 //last age
  init_int compfage                             //first age for composition data
  init_int complage                             //last age for composition data

  //Fishery 
  init_vector obs_catch_flt1(fyear,lyear)       //Catch Fleet 1
  init_number SE_catch_flt1                     //SE of catch for each year in fleet 1

  //Survey
  init_int nyr_Acoustic                             //number of years with CPUE data, same as fleet 2
  init_vector Acoustic_yrs(1,nyr_Acoustic)              //Indices of years with CPUE data, same as fleet 2
  init_vector Acoustic_data(1,nyr_Acoustic)             //CPUE data 
  init_vector Acoustic_SE(1,nyr_Acoustic)               //CPUE data SEs for each yr

  //Composition Data
  //Fishery (Fleet 1)
  init_int fyr_flt1                             //first year fleet 1
  init_int lyr_flt1                             //last year fleet 1
  init_int month_flt1                           //Month of fleet 1
  init_int AE_defstart_flt1                     //Start index for aging error fleet 1
  init_int AE_defend_flt1                     //end index for aging error fleet 1
  init_vector SS_flt1(fyr_flt1,lyr_flt1)              //sample size for fleet 1
  init_matrix ObsComp_flt1(fyr_flt1,lyr_flt1,compfage,complage)  //observed comp flt 1, contains non-proportions, corrected later

  //Survey (Fleet 2)
  init_int month_flt2                           //Month of fleet 2
  init_vector AE_defyrs_flt2(1,nyr_Acoustic)        //Year indices for aging error fleet 2
  init_vector SS_flt2(1,nyr_Acoustic)               //sample size for fleet 2
  init_matrix ObsComp_flt2(1,nyr_Acoustic,compfage,complage)  //observed comp flt 2, , contains non-proportions, corrected later

  //Aging Error
  init_int num_AE_def                           //number of aging error definitions
  init_matrix AE_means(1,num_AE_def,fage,lage)            //Aging error means
  init_matrix AE_SEs(1,num_AE_def,fage,lage)              //Aging error SEs
  !!AE_means=AE_means-0.5;         //Fixing AE_means, i think

  //Wt-at-Age
  //Weight at age at the start of the year
  init_int WAA1_preyear_index
  init_vector WAA1_preyear(fage,lage)
  init_int WAA1_fyear_index
  init_int WAA1_lyear_index
  init_matrix WAA1(WAA1_fyear_index,WAA1_lyear_index,fage,lage)

  //Weight at age in the middle of the year
  init_int WAA7_preyear_index
  init_vector WAA7_preyear(fage,lage)
  init_int WAA7_fyear_index
  init_int WAA7_lyear_index
  init_matrix WAA7(WAA7_fyear_index,WAA7_lyear_index,fage,lage)

  //Weight at age for Fishery (fleet 1)
  init_int WAAflt1_preyear_index
  init_vector WAAflt1_preyear(fage,lage)
  init_int WAAflt1_fyear_index
  init_int WAAflt1_lyear_index
  init_matrix WAAflt1(WAAflt1_fyear_index,WAAflt1_lyear_index,fage,lage)

  //Weight at age for Survey (fleet 2)
  init_int WAAflt2_preyear_index
  init_vector WAAflt2_preyear(fage,lage)
  init_int WAAflt2_fyear_index
  init_int WAAflt2_lyear_index
  init_matrix WAAflt2(WAAflt2_fyear_index,WAAflt2_lyear_index,fage,lage)
  
  //Maturity*Fecundity
  init_int MatFec_preyear_index
  init_vector MatFec_preyear(fage,lage)
  init_int MatFec_fyear_index
  init_int MatFec_lyear_index
  init_matrix MatFec(MatFec_fyear_index,MatFec_lyear_index,fage,lage)

  init_int test
  int i
  int j
  int k
  int mycount
  !!mycount=0;
  
//  !!cout << test << endl;
//  !!exit(88);
 
  vector ones(compfage,complage)                         //ones for DM     
  !!ones=1;
  vector fsel(fage,lage)
  matrix mat_fsel(fyear,lyear,fage,lage)
  matrix obs_fcomp(fyr_flt1,lyr_flt1,compfage,complage)  //observed comp flt 1
  matrix obs_survcomp(1,nyr_Acoustic,compfage,complage)  //observed comp flt 2

  int recdev_start
  !!recdev_start=1946;
  int recdev_end
  !!recdev_end=2020;
  vector b_recbias(recdev_start,recdev_end) //Taylor and Methot recbias ramp
  number bmax                                           //maximum value for bias ramp
  !!bmax=0.87;
  number y1_b
  !!y1_b=1965;                                          //first year of the bias ramp up adjustment period
  number y2_b
  !!y2_b=1971;                                          //last year of the bias ramp up adjustment period
  number y3_b
  !!y3_b=2018;                                          //first year of the bias ramp down adjustment period
  number y4_b
  !!y4_b=2019;                                          //last year of the bias ramp down adjustment period

INITIALIZATION_SECTION
  //   Have it all in pin file
  
PARAMETER_SECTION
  //Population
  init_bounded_number log_R0(13.,17.,1)                            //Unfished recruitment
  init_bounded_number steepness(0.2,1,4)                           //Steepness SR curve
  init_bounded_vector log_rec_devs(fyear-lage,recdev_end,-6.,6.,2) //recruitment deviations
  init_bounded_number log_M(-3,-0.9,4)                             //Natural Mortality
  //Fishery
  init_bounded_vector log_Fint(fyear,lyear,-10.0,0.0,1)            //Fully selected fishing mortalities (will code in hybrid method later
  init_bounded_vector fsel_par(2,6,-5.,9.,2)                       //Non-parametric fishery selectivity
  init_bounded_vector sel_devs(1,150,-5.,5.,5)                     //Selectivity devs, matrix from years 1991-2020 from ages 2-6
  init_bounded_number log_ftheta(-5.,20.,5)                        //DM overdispersion parameter
  //Survey
  init_bounded_vector survsel_par(3,6,-5.,9.,2)                    //Non-parametric survey selectivity 
  init_bounded_number log_surv_sd(-3.,0.2,5)                       //Additional survey sd
  init_bounded_number log_survtheta(-5.,20.,5)                     //DM overdispersion parameter

 //Derived variables
  matrix N(fyear,lyear+1,fage,lage)
  matrix N7(fyear,lyear+1,fage,lage)                               //Abundance at age in the middle of the year
  matrix F(fyear,lyear,fage,lage)
  matrix M(fyear,lyear,fage,lage)
  matrix Z(fyear,lyear,fage,lage)

 //Sel Vars
  vector pa_fishery(fage,lage)
  matrix pa_fishery_time(fyear,lyear,fage,lage)
  matrix Saprime_fishery(fyear,lyear,fage,lage)
  matrix Sel_fishery(fyear,lyear,fage,lage)
  matrix sel_devs_mat(1991,2020,2,6)

  vector pa_surv(fage,lage)
  vector pa_surv_time(fage,lage)
  vector Saprime_surv(fage,lage)
  vector Sel_surv(fage,lage)

  vector lxo(fage,lage)                     //unfished survivorship
  vector N0_age(fage,lage)                      //Unfished abundance at age
  number SSB0                               //Numbers at age * Fecundity
  vector SpawnBio(fyear,lyear+1)              //Spawning biomass

  //Aging Error
  3darray ageing_error(1,num_AE_def,fage,lage,fage,lage)       //aging error matrix
  3darray ageerr_crop(1,num_AE_def,fage,lage,compfage,complage)
  
  //Fishery 
  matrix preErr_pred_caa(fyear,lyear,fage,lage) //Predicted Catch at age
  matrix pred_caa(fyear,lyear,fage,lage)        //Predicted Catch at age with aging error
  vector pred_harv(fyear,lyear)                 //Predicted Harvest 
  matrix pred_fcomp(fyear,lyear,fage,lage)      //Predicted Fishery Composition
  matrix pred_fcomp_compressed(fyear,lyear,compfage,complage)      //Predicted Fishery Composition

  //Survey
  vector Acoustic_sd(1,nyr_Acoustic)
  matrix preErr_pred_Survaa(fyear,lyear,fage,lage)
  matrix pred_Survaa(fyear,lyear,fage,lage)
  vector pred_Survindex(fyear,lyear)
  matrix pred_Survcomp(fyear,lyear,fage,lage)
  matrix pred_Survcomp_compressed(fyear,lyear,compfage,complage)
  vector pred_SurvBio(fyear,lyear)
  vector pred_SurvBio_ind(1,nyr_Acoustic)
  sdreport_number log_survq                               //Catchability for acoustic survey (analytical estimate as median unbiased)

  number ftheta                                  //DM overdispersion parameter
  vector fBeta(fyr_flt1,lyr_flt1)                //DM overdispersion parameter
  number survtheta                               //DM overdispersion parameter
  vector survBeta(1,nyr_Acoustic)                //DM overdispersion parameter

  number pi
  number ftheta_sig                              //Sigma for theta prior
  number survtheta_sig                           //Sigma for theta prior
  number M_sd                                    //Variance of Natural mortality prior
  number log_priorM                              //Prior for M
  number alpha_h                                 //Steepness prior alpha (for Beta dist)
  number beta_h                                  //Steepness prior beta  (for Beta dist)
  number seldevs_sd                              //sd for seldevs prior
  number log_rec_sd                              //sd for log recruitment

  number mean_beta
  number sd_beta
  number LB_beta
  number UB_beta
  number v_beta
  number tau_beta
  number m_beta
  number n_beta

  vector exp_recdevs(recdev_start,recdev_end+1) //exponentiated rec_devs and corrected for bias 

//Likelihood components
  number NLLcatch
  number NLL1
  number NLL2
  number NLL3
  number NLL4
  number NLL5
  number NLL6
  number NLL7
  number NLL8
  number NLL9

  objective_function_value NLL;

PRELIMINARY_CALCS_SECTION

//Make Composition data bc it's currently not in prop form
  for(i=fyr_flt1;i<=lyr_flt1;i++){
   obs_fcomp(i)=ObsComp_flt1(i)/sum(ObsComp_flt1(i));
   obs_fcomp(i)=(obs_fcomp(i)+0.001)/sum(obs_fcomp(i)+0.001);
  }

  for(j=1;j<=nyr_Acoustic;j++){
   obs_survcomp(j)=ObsComp_flt2(j)/sum(ObsComp_flt2(j));
   obs_survcomp(j)=(obs_survcomp(j)+0.001)/sum(obs_survcomp(j)+0.001);
  }
  
  for(i=0;i<=29;i++){   //1991-2020
   for(j=1;j<=5;j++){   //ages 3-7
    sel_devs_mat(i+1991,j+1)=sel_devs(i*5+j);
   }
  }

PROCEDURE_SECTION
  get_ageing_error();
  get_selectivity();
  get_mortality();
  get_population();
  get_catch();
  get_acoustic();
  get_objective();
  
  if (mceval_phase())
    {
    ofstream mcres("Params_Hake.txt",ios::app);
    if (mycount==0){
     mcres << "log_R0 steepness log_M log_surv_sd log_ftheta log_survtheta \
             fsel_par2 fsel_par3 fsel_par4 fsel_par5 fsel_par6 \
             survsel_par3 survsel_par4 survsel_par5 survsel_par6 \
             sel_dev2_91 sel_dev3_91 sel_dev4_91 sel_dev5_91 sel_dev6_91 \
             sel_dev2_92 sel_dev3_92 sel_dev4_92 sel_dev5_92 sel_dev6_92 \
             sel_dev2_93 sel_dev3_93 sel_dev4_93 sel_dev5_93 sel_dev6_93 \
             sel_dev2_94 sel_dev3_94 sel_dev4_94 sel_dev5_94 sel_dev6_94 \
             sel_dev2_95 sel_dev3_95 sel_dev4_95 sel_dev5_95 sel_dev6_95 \
             sel_dev2_96 sel_dev3_96 sel_dev4_96 sel_dev5_96 sel_dev6_96 \
             sel_dev2_97 sel_dev3_97 sel_dev4_97 sel_dev5_97 sel_dev6_97 \
             sel_dev2_98 sel_dev3_98 sel_dev4_98 sel_dev5_98 sel_dev6_98 \
             sel_dev2_99 sel_dev3_99 sel_dev4_99 sel_dev5_99 sel_dev6_99 \
             sel_dev2_00 sel_dev3_00 sel_dev4_00 sel_dev5_00 sel_dev6_00 \
             sel_dev2_01 sel_dev3_01 sel_dev4_01 sel_dev5_01 sel_dev6_01 \
             sel_dev2_02 sel_dev3_02 sel_dev4_02 sel_dev5_02 sel_dev6_02 \
             sel_dev2_03 sel_dev3_03 sel_dev4_03 sel_dev5_03 sel_dev6_03 \
             sel_dev2_04 sel_dev3_04 sel_dev4_04 sel_dev5_04 sel_dev6_04 \
             sel_dev2_05 sel_dev3_05 sel_dev4_05 sel_dev5_05 sel_dev6_05 \
             sel_dev2_06 sel_dev3_06 sel_dev4_06 sel_dev5_06 sel_dev6_06 \
             sel_dev2_07 sel_dev3_07 sel_dev4_07 sel_dev5_07 sel_dev6_07 \
             sel_dev2_08 sel_dev3_08 sel_dev4_08 sel_dev5_08 sel_dev6_08 \
             sel_dev2_09 sel_dev3_09 sel_dev4_09 sel_dev5_09 sel_dev6_09 \
             sel_dev2_10 sel_dev3_10 sel_dev4_10 sel_dev5_10 sel_dev6_10 \
             sel_dev2_11 sel_dev3_11 sel_dev4_11 sel_dev5_11 sel_dev6_11 \
             sel_dev2_12 sel_dev3_12 sel_dev4_12 sel_dev5_12 sel_dev6_12 \
             sel_dev2_13 sel_dev3_13 sel_dev4_13 sel_dev5_13 sel_dev6_13 \
             sel_dev2_14 sel_dev3_14 sel_dev4_14 sel_dev5_14 sel_dev6_14 \
             sel_dev2_15 sel_dev3_15 sel_dev4_15 sel_dev5_15 sel_dev6_15 \
             sel_dev2_16 sel_dev3_16 sel_dev4_16 sel_dev5_16 sel_dev6_16 \
             sel_dev2_17 sel_dev3_17 sel_dev4_17 sel_dev5_17 sel_dev6_17 \
             sel_dev2_18 sel_dev3_18 sel_dev4_18 sel_dev5_18 sel_dev6_18 \
             sel_dev2_19 sel_dev3_19 sel_dev4_19 sel_dev5_19 sel_dev6_19 \
             sel_dev2_20 sel_dev3_20 sel_dev4_20 sel_dev5_20 sel_dev6_20 \
             log_rec_dev46 log_rec_dev47 log_rec_dev48 log_rec_dev49 log_rec_dev50 log_rec_dev51 log_rec_dev52 log_rec_dev53 log_rec_dev54 log_rec_dev55 log_rec_dev56 log_rec_dev57 log_rec_dev58 log_rec_dev59 log_rec_dev60 log_rec_dev61 log_rec_dev62 log_rec_dev63 log_rec_dev64 log_rec_dev65 log_rec_dev66 log_rec_dev67 log_rec_dev68 log_rec_dev69 log_rec_dev70 log_rec_dev71 log_rec_dev72 log_rec_dev73 log_rec_dev74 log_rec_dev75 log_rec_dev76 log_rec_dev77 log_rec_dev78 log_rec_dev79 log_rec_dev80 log_rec_dev81 log_rec_dev82 log_rec_dev83 log_rec_dev84 log_rec_dev85 log_rec_dev86 log_rec_dev87 log_rec_dev88 log_rec_dev89 log_rec_dev90 log_rec_dev91 log_rec_dev92 log_rec_dev93 log_rec_dev94 log_rec_dev95 log_rec_dev96 log_rec_dev97 log_rec_dev98 log_rec_dev99 log_rec_dev00 log_rec_dev01 log_rec_dev02 log_rec_dev03 log_rec_dev04 log_rec_dev05 log_rec_dev06 log_rec_dev07 log_rec_dev08 log_rec_dev09 log_rec_dev10 log_rec_dev11 log_rec_dev12 log_rec_dev13 log_rec_dev14 log_rec_dev15 log_rec_dev16 log_rec_dev17 log_rec_dev18 log_rec_dev19 log_rec_dev20 \
             log_Fint_66 log_Fint_67 log_Fint_68 log_Fint_69 log_Fint_70 log_Fint_71 log_Fint_72 log_Fint_73 log_Fint_74 log_Fint_75 log_Fint_76 log_Fint_77 log_Fint_78 log_Fint_79 log_Fint_80 log_Fint_81 log_Fint_82 log_Fint_83 log_Fint_84 log_Fint_85 log_Fint_86 log_Fint_87 log_Fint_88 log_Fint_89 log_Fint_90 log_Fint_91 log_Fint_92 log_Fint_93 log_Fint_94 log_Fint_95 log_Fint_96 log_Fint_97 log_Fint_98 log_Fint_99 log_Fint_00 log_Fint_01 log_Fint_02 log_Fint_03 log_Fint_04 log_Fint_05 log_Fint_06 log_Fint_07 log_Fint_08 log_Fint_09 log_Fint_10 log_Fint_11 log_Fint_12 log_Fint_13 log_Fint_14 log_Fint_15 log_Fint_16 log_Fint_17 log_Fint_18 log_Fint_19 log_Fint_20"<<endl;
     }
     mcres << log_R0 << " " << steepness << " " << log_M << " " << log_surv_sd << " " << log_ftheta << " " << log_survtheta << " "
           << fsel_par << " "
           << survsel_par << " "
           << sel_devs << " "
           << log_rec_devs << " "
           << log_Fint <<  endl;
     mycount++;
    }
    
FUNCTION get_ageing_error
//Ageing Error
 //If there is no error, identity matrix
  if (age_err==0){
   for(k=1;k<=num_AE_def;k++){    //k represents year
    ageing_error(k)=identity_matrix(fage,lage);
   }
  }else if (age_err==1){
//SS uses the normal distribution, so ageing error matrix every year represents the probability that true age j gets coded age i
   for(k=1;k<=num_AE_def;k++){    //k represents year
    for(j=fage;j<=lage;j++){     //j represents true age
     ageing_error(k,j,fage)=cumd_norm(((int(fage)+0.5)-AE_means(k,j))/AE_SEs(k,j));
     for(i=fage+1;i<=lage-1;i++){  //i represents coded age 
      ageing_error(k,j,i)=cumd_norm(((int(i)+0.5)-AE_means(k,j))/AE_SEs(k,j)) - cumd_norm(((int(i)-0.5)-AE_means(k,j))/AE_SEs(k,j));  
     }
     ageing_error(k,j,lage)=1-cumd_norm(((int(lage)-0.5)-AE_means(k,j))/AE_SEs(k,j));
    }
   }
  }
   //Making first age to be classified as age 1, so no age 1s are classified as zero
  for(k=1;k<=num_AE_def;k++){    //k represents year
   for(j=fage;j<=lage;j++){     //j represents true age
    ageing_error(k,j,1)+=ageing_error(k,j,0);
    ageing_error(k,j,0)=0; //making coded age zero column zero
   }
  }

//This is just code that compresses the age err matrix to make it start at age 1 to plus group
  /*
  for(k=1;k<=num_AE_def;k++){    //k represents year
   for(j=fage;j<=lage;j++){     //j represents true age
    ageerr_crop(k,j,1)=sum(ageing_error(k,j)(0,1));
    ageerr_crop(k,j)(2,14)=ageing_error(k,j)(2,14);
    ageerr_crop(k,j,15)=sum(ageing_error(k,j)(15,20));
   }
  }
  */

FUNCTION get_selectivity
//Fishery Selectivity
   //Non-parametric ages 1-6, then fixed at 6 after that (figure 22 in assessment report)
   //And time varying
  pa_fishery(0)=-1000;
  pa_fishery(1)=0;
  pa_fishery(2,6)=fsel_par;
  pa_fishery(7,lage)=0;

  for(i=fyear;i<=lyear;i++){
   for(j=fage;j<=lage;j++){
    pa_fishery_time(i,j) = pa_fishery(j);
     if(i>=1991 & j>=2 & j<7){pa_fishery_time(i,j)=pa_fishery(j)+sel_devs_mat(i,j);}
   }
  }
    
  for(i=fyear;i<=lyear;i++){
   for(j=fage+1;j<=lage;j++){
    Saprime_fishery(i,j)=sum(pa_fishery_time(i)(fage+1,j));
   }
  }
  
  for(i=fyear;i<=lyear;i++){
   for(j=fage+1;j<=lage;j++){
    Sel_fishery(i,j)=mfexp(Saprime_fishery(i,j)-max(Saprime_fishery(i)(fage+1,lage)));
   }
  }
  
  //Alternative Time-varying Selectivity, with exponentiation
  /*
  for(i=fyear;i<=lyear;i++){
   for(j=fage;j<=lage;j++){
    pa_fishery_time(i,j) = pa_fishery(j);
   }
  }
  for(i=fyear;i<=lyear;i++){
    for(j=fage;j<=lage;j++){
    Saprime_fishery(i,j)=sum(pa_fishery_time(i)(fage,j));
   }
  }
  for(i=fyear;i<=lyear;i++){
   Sel_fishery(i,0)=0.0;
    for(j=fage+1;j<=lage;j++){
     Sel_fishery(i,j)=mfexp(Saprime_fishery(i,j)-max(Saprime_fishery(i)));
     if (i>=1991 & j>=2 & j<7){ Sel_fishery(i,j)=mfexp(Saprime_fishery(i,j)-max(Saprime_fishery(i)))*mfexp(sel_devs_mat(i,j));}
    }
  }
  */

//Acoustic Survey Selectivity
  //Ages 2-5, fixed after that at age 6 level (Figure 24 assessment report)
  pa_surv(0,1)=-1000;
  pa_surv(2)=0;
  pa_surv(3,6)=survsel_par;
  pa_surv(7,lage)=0;

  for(j=fage+2;j<=lage;j++){
   Saprime_surv(j)=sum(pa_surv(fage+2,j));
  } 

  for(j=fage+2;j<=lage;j++){
   Sel_surv(j)=mfexp(Saprime_surv(j)-max(Saprime_surv(2,lage)));
  }

//  cout<<Sel_surv<<endl;
//  cout<<" "<<endl;
//  cout<<Sel_fishery<<endl;
//  exit(88);
  
FUNCTION get_mortality
  //MORTALITY      //Will code in Hybrid method later, currently Fint as params
  for(i=fyear;i<=lyear;i++){
   for(j=fage;j<=lage;j++){
    F(i,j)=Sel_fishery(i,j)*mfexp(log_Fint(i));  
    M(i,j)=mfexp(log_M);         
   }
  }
  Z=M+F;         //Total instantaneous mortality

//  cout << Z << endl;
//  exit(88);
  
FUNCTION get_population
 //Getting bias ramp, and the exponentiated rec devs
  if (b_ramp==0){
    b_recbias=1;
   for(i=recdev_start;i<=recdev_end;i++){
    exp_recdevs(i)=mfexp(log_rec_devs(i)-0.5*square(mfexp(log_rec_sd)));
   }
  } else if (b_ramp==1){
  for(i=recdev_start;i<=recdev_end;i++){
   if(i>y1_b & i<y2_b){        b_recbias(i) = bmax*((i-y1_b)/(y2_b-y1_b));}
   else if(i>=y2_b & i<=y3_b){ b_recbias(i) = bmax;}
   else if(i>y3_b & i<y4_b) {  b_recbias(i) = bmax*(1-(y3_b-i)/(y4_b-y3_b));}
   else if(i<=y1_b | i>=y4_b){ b_recbias(i) = 0;}
    
   exp_recdevs(i)=mfexp(-0.5*b_recbias(i)*square(mfexp(log_rec_sd))+log_rec_devs(i));
  }}
   exp_recdevs(recdev_end+1)=1;               // Recdev for 2021 recruitment just mean rec. 

 //Unfished Spawning biomass calcs
  lxo(fage)=1;
  for(j=fage+1;j<=lage;j++){
   lxo(j)=lxo(j-1)*mfexp(-1.0*mfexp(log_M));   //cumulative survival
  }
  lxo(lage)=lxo(lage)/(1-mfexp(-1.0*mfexp(log_M)));
  N0_age=mfexp(log_R0)*lxo;
  SSB0=sum(elem_prod(N0_age,MatFec_preyear));   //Numbers at age * Fecundity

  N.initialize();
  //Unfished Abundance at age
  for(j=fage;j<=lage;j++){
   N(fyear,j)=mfexp(log_R0)*exp_recdevs(1966-j)*lxo(j);
  }
  //Spawning Biomass prior the first year
  SpawnBio(fyear)=sum(elem_prod(N(fyear),MatFec_preyear));  // Getting spawning biomass for the first year (acounting for rec dev, unlike SSB0_FLA)

 //Population loop for model time series
  for(i=fyear+1;i<=lyear+1;i++){
    for(j=fage+1;j<=lage;j++)   {
     N(i,j)=N(i-1,j-1)*mfexp(-1.0*Z(i-1,j-1));
    }
   N(i,lage)+=N(i-1,lage)*mfexp(-1.0*Z(i-1,lage));   //Plus group
    if (i<MatFec_fyear_index){ SpawnBio(i)=sum(elem_prod(N(i),MatFec_preyear)); } //Spawning Biomass
     else if(i>=MatFec_fyear_index){ SpawnBio(i)=sum(elem_prod(N(i),MatFec(i))); }
    N(i,fage)=((4.*steepness*mfexp(log_R0)*SpawnBio(i))/(SSB0*(1.-steepness)+SpawnBio(i)*(5.*steepness-1.)))*exp_recdevs(i);         //Recruitment with bias correction
  }

//  cout<< lxo <<endl;
//  exit(1);

FUNCTION get_catch
  ///*
//Fishery, Fleet 1
  for (i=fyear;i<=lyear;i++){ 
   preErr_pred_caa(i)=elem_prod(elem_div(F(i),Z(i)),elem_prod(1.0-mfexp(-1.0*Z(i)),N(i)));      //Baranov catch equation for predicting catch at age
   if(i>=fyr_flt1){pred_caa(i)=preErr_pred_caa(i)*ageing_error(i-(fyr_flt1-AE_defstart_flt1));} //Applying ageing error to catch at age
   else if(i<fyr_flt1){pred_caa(i)=preErr_pred_caa(i);}                                         //Applying ageing error to catch at age
   pred_fcomp(i)=pred_caa(i)/sum(pred_caa(i));                                                  //calculating predicted catch age composition
   pred_fcomp_compressed(i)(compfage,complage-1)=pred_fcomp(i)(compfage,complage-1);            //correcting for observed ages 
   pred_fcomp_compressed(i,complage)=sum(pred_fcomp(i)(complage,lage));                         //correcting for observed ages (plus group)
   pred_fcomp_compressed(i)=(pred_fcomp_compressed(i)+0.001)/sum(pred_fcomp_compressed(i)+0.001);                                  //adding small constant
   if(i>=fyr_flt1){pred_harv(i)=sum(elem_prod(preErr_pred_caa(i),WAAflt1(i)));}                                  //Predicted total harvest weight each year
   else if(i<fyr_flt1){pred_harv(i)=sum(elem_prod(preErr_pred_caa(i),WAAflt1_preyear));}                                  //Predicted total harvest weight each year
  }
  //*/
  
//  cout << pred_fcomp << endl;
//  cout << " " << endl;
//  cout << pred_fcomp_compressed << endl;
//  exit(88);

FUNCTION get_acoustic
  ///*
   //Acoustic Survey Index and Composition
  for (i=fyear;i<=lyear;i++){
   N7(i) = elem_prod(N(i),mfexp(-0.5*Z(i)));                                       //abundance at age in the middle of the year
   preErr_pred_Survaa(i)=elem_prod(Sel_surv,N7(i));	                                  //Predicted survey caa (before aging error), correct for season
   if(i<1975){pred_SurvBio(i)=sum(elem_prod(preErr_pred_Survaa(i),WAA7_preyear));}            //Predicted biomass to use for survey index
   else if(i>=1975){pred_SurvBio(i)=sum(elem_prod(preErr_pred_Survaa(i),WAA7(i)));}       //Predicted biomass to use for survey index
   if(i<1973){pred_Survaa(i)=preErr_pred_Survaa(i);}                                      //No ageing error prior to 1973
   else if(i>=1973){pred_Survaa(i)=preErr_pred_Survaa(i)*ageing_error(i-1972);}           //Applying ageing error to survey caa
   pred_Survcomp(i)=pred_Survaa(i)/sum(pred_Survaa(i));                                   //Predicted Survey composition
   pred_Survcomp_compressed(i)(compfage,complage-1)=pred_Survcomp(i)(compfage,complage-1);  //correcting for observed ages 
   pred_Survcomp_compressed(i,complage)=sum(pred_Survcomp(i)(complage,lage));               //correcting for observed ages (plus group)
   pred_Survcomp_compressed(i)=(pred_Survcomp_compressed(i)+0.001)/sum(pred_Survcomp_compressed(i)+0.001);               //adding small constant
  }
  //*/

//Subsetting for only the years which have index data so I can calculate q
  for(i=fyear;i<=lyear;i++){
   for(j=1;j<=nyr_Acoustic;j++){
    if(i==Acoustic_yrs(j)){
     pred_SurvBio_ind(j)=pred_SurvBio(i);                                                 //Subsetting for only the years which have index data so I can calculate q
    }
   }
  }

//Analytical log catchability as the exp( mean(log(obs_index)) )
  log_survq = mean(log(elem_div(Acoustic_data,pred_SurvBio_ind)));                        //Getting analytical q
  pred_Survindex=mfexp(log_survq)*pred_SurvBio;                                    //Predicted survey index

//  cout << pred_Survcomp << endl;
//  cout << " " << endl;
//  cout << pred_Survcomp_compressed << endl;
//  exit(88);
 
FUNCTION get_objective

////////////////////////////////////
//LIKELIHOODS for data components
////////////////////////////////////
  pi=3.141593;

//Likelihood for Fishery Catch, Lognormal, Needed because I am estimating F as a parameter
  //NLLcatch = sum(log(obs_catch_flt1))+double(size_count(obs_catch_flt1))*0.5*log(2*pi)+double(size_count(obs_catch_flt1))*log(SE_catch_flt1)+norm2(log(obs_catch_flt1)-log(pred_harv))/(2*square(SE_catch_flt1));
  NLLcatch = norm2(log(obs_catch_flt1)-log(pred_harv+1.0E-6))/(2*square(SE_catch_flt1));

//Likelihood for Fishery Composition, Dirichlet Multinomial 
  ftheta=mfexp(log_ftheta);
  fBeta=ftheta*SS_flt1;
  NLL1.initialize();
  for(i=fyr_flt1;i<=lyr_flt1;i++){
  //Linear Version from Thorson
   if(SS_flt1(i)>0){
    NLL1 += -1*(gammln(SS_flt1(i)+1)-sum(gammln(SS_flt1(i)*obs_fcomp(i)+ones))+gammln(ftheta)-gammln(SS_flt1(i)+ftheta)+sum(gammln(SS_flt1(i)*obs_fcomp(i)+ftheta*pred_fcomp_compressed(i))-gammln(ftheta*pred_fcomp_compressed(i)))); 
   }
  }
  
//Likelihood for Acoustic Survey Index, lognormal likelihood
  Acoustic_sd=Acoustic_SE+mfexp(log_surv_sd);    //for additive sd's
//  Acoustic_sd=sqrt(square(Acoustic_SE)+square(mfexp(log_surv_sd)));    //for additive variances
  NLL2.initialize();
  for(i=fyear;i<=lyear;i++){
   for(j=1;j<=nyr_Acoustic;j++){
    if(i==Acoustic_yrs(j)){
     //NLL2 += log(Acoustic_data(j))+0.5*log(2*pi)+log(Acoustic_sd(j))+square(log(Acoustic_data(j))-log(pred_Survindex(i)))/(2*square(Acoustic_sd(j)));
     NLL2 += log(Acoustic_sd(j))+square(log(Acoustic_data(j))-log(pred_Survindex(i)))/(2*square(Acoustic_sd(j)));
     //NLL2 += Acoustic_SE(j)*log_surv_sd+square(log(Acoustic_data(j))-log(pred_Survindex(i)))/(2*square(Acoustic_sd(j)));
    }
   }
  }

//Likelihood for Survey Composition, Dirichlet-Multinomial, from ages 2-15
  survtheta=mfexp(log_survtheta);
  survBeta=survtheta*SS_flt2;
  NLL3.initialize();
  for(i=fyear;i<=lyear;i++){
   for(j=1;j<=nyr_Acoustic;j++){
    if(i==Acoustic_yrs(j)){
     //Linear Version from Thorson
     if(SS_flt2(j)>0){
      NLL3 += -1*(gammln(SS_flt2(j)+1)-sum(gammln(SS_flt2(j)*obs_survcomp(j)(compfage+1,complage)+ones(compfage+1,complage)))+gammln(survtheta)-gammln(SS_flt2(j)+survtheta)+sum(gammln(SS_flt2(j)*obs_survcomp(j)(compfage+1,complage)+survtheta*pred_Survcomp_compressed(i)(compfage+1,complage))-gammln(survtheta*pred_Survcomp_compressed(i)(compfage+1,complage))));
     }
    }
   }
  }

////////////////////////
//Priors
////////////////////////

  //Prior on Natural Mortality, Lognormal(0.20,1.11)
  M_sd=0.1;
  log_priorM=log(0.20);
  //NLL4 = log_M+0.5*log(2*pi)+log(M_sd)+square(log_M-log_priorM)/(2*square(M_sd));
  NLL4 = square(log_M-log_priorM)/(2*square(M_sd));

  //This is the actual Beta Negative Log-Likelihood, however SS uses formulation further below
  /*
  //Prior on Steepness, Beta(9.76,2.80), this equals mean=0.78 and sd=0.11
  alpha_h=9.76;      
  beta_h=2.80;
  NLL5 = (alpha_h-1)*steepness+(beta_h+1)*log(1-steepness)-log((mfexp(gammln(alpha_h))*mfexp(gammln(beta_h)))/mfexp(gammln(alpha_h+beta_h)));
  */
  
  //Beta Negative log-likelihood from CASAL and from Stock synthesis, they are the same (on a relative scale)
  mean_beta=0.777;                                                    //Prior mean for Beta
  sd_beta=0.113;                                                      //Prior sd for Beta
  LB_beta=0.2;                                                        //Lower bound for Beta
  UB_beta=1;                                                          //Upper bound for Beta
  v_beta=(mean_beta-LB_beta)/(UB_beta-LB_beta);
  tau_beta=(((mean_beta-LB_beta)*(UB_beta-mean_beta))/square(sd_beta))-1;
  m_beta=tau_beta*v_beta;
  n_beta=tau_beta*(1-v_beta);
  //NLL5 = (1-m_beta)*log(steepness-LB_beta)+(1-n_beta)*log(UB_beta-steepness);  //Beta from CASAL manual
  NLL5 = (1-m_beta)*log(steepness-LB_beta)+(1-n_beta)*log(UB_beta-steepness)-(1-m_beta)*log(mean_beta-LB_beta)-(1-n_beta)*log(UB_beta-mean_beta);  //Beta from SS Manual
 
  //Recruitment deviations prior Lognormal(0,rec_sigma)
  log_rec_sd=log(1.4);
  //NLL6 = sum(log_rec_devs)+double(size_count(log_rec_devs))*0.5*log(2*pi)+0.5*sum((norm2(log_rec_devs)/square(mfexp(log_rec_sd)))+b_recbias*log_rec_sd); 
  NLL6 = 0.5*sum((square(log_rec_devs)/square(mfexp(log_rec_sd)))+b_recbias*log(square(mfexp(log_rec_sd)))); 

  //Selectivity deviations prior Normal(0,1.4)
  seldevs_sd=1.4;
  //NLL7 = double(size_count(sel_devs))*0.5*log(2*pi)+double(size_count(sel_devs))*log(seldevs_sd)+norm2(sel_devs)/(2*square(seldevs_sd));
  NLL7 = norm2(sel_devs)/(2*square(seldevs_sd));

  //DM Overdispersion param Fishery, Normal(0,1.813)
  ftheta_sig=1.813;
  //NLL8 = 0.5*log(2*pi)+log(ftheta_sig)+square(log_ftheta)/(2*square(ftheta_sig));
  NLL8 = square(log_ftheta)/(2*square(ftheta_sig));

  //DM Overdispersion param Survey, Normal(0,1.813)
  survtheta_sig=1.813;
  //NLL9 = 0.5*log(2*pi)+log(survtheta_sig)+square(log_survtheta)/(2*square(survtheta_sig));
  NLL9 = square(log_survtheta)/(2*square(survtheta_sig));

//Total NLL
  NLL=NLL1+NLL2+NLL3+NLL4+NLL5+NLL6+NLL7+NLL8+NLL9+NLLcatch;

  /*
  cout << NLL << endl;
  cout << NLLcatch << endl;
  cout << NLL1 << endl;
  cout << NLL2 << endl;
  cout << NLL3 << endl;
  cout << NLL4 << endl;
  cout << NLL5 << endl;
  cout << NLL6 << endl;
  cout << NLL7 << endl;
  cout << NLL8 << endl;
  cout << NLL9 << endl;
  exit(1);
  */  

GLOBALS_SECTION
  #include "admodel.h"
  #include "admb2r.cpp"

REPORT_SECTION
  report<<"Final gradient: "<<objective_function_value::pobjfun->gmax << endl;
  report << " " << endl;
  report <<"NLL   "<< NLL<<endl;
  report << " " << endl;
  report << "Obs Harvest" << "   " << "Pred Harvest" << endl;
  for(i=fyear;i<=lyear;i++){
    report << obs_catch_flt1(i) << "         " << pred_harv(i) << endl;
   }
  report << " " << endl;
  report << "Obs Acoustic" << "   " << "Pred Acoustic" << endl;
  for(i=fyear;i<=lyear;i++){
   for(j=1;j<=nyr_Acoustic;j++){
    if(i==Acoustic_yrs(j)){
    report << Acoustic_data(j) << "        " << pred_Survindex(i) << endl;
   }}}
  report << " " << endl;
  report << "Obs Surv Comp " << endl;
  report << obs_survcomp << endl;
  report << " " << endl;
  report << "Pred Surv Comp " << endl;
  report << pred_Survcomp_compressed << endl;
  report << " " << endl;
  report << "Obs Fishery Comp " << endl;
  report << obs_fcomp << endl;
  report << " " << endl;
  report << "Pred Fishery Comp " << endl;
  report << pred_fcomp_compressed << endl;
  report << " " << endl;
  report << "Fishery Selectivity" << endl;
  report << Sel_fishery << endl;
  report << " " << endl;
  report << "Survey Selectivity" << endl;
  report << Sel_surv << endl;
  report << " " << endl;
  report << NLLcatch << endl;
  report << NLL1 << endl;
  report << NLL2 << endl;
  report << NLL3 << endl;
  report << NLL4 << endl;
  report << NLL5 << endl;
  report << NLL6 << endl;
  report << NLL7 << endl;
  report << NLL8 << endl;
  report << NLL9 << endl;

//Writing results into R
  open_r_file("PHake.rdat");
   open_r_vector("Gradient_NLL");
    wrt_r_item("Gradient",objective_function_value::pobjfun->gmax);
    wrt_r_item("NLL",NLL);
   close_r_vector();
    open_r_matrix("N");
      wrt_r_matrix(N);
    close_r_matrix();
    wrt_r_complete_vector("N0_age",N0_age);
    wrt_r_complete_vector("SpawnBio",SpawnBio);
    wrt_r_complete_vector("Predharv", pred_harv);
    wrt_r_complete_vector("Obsharv", obs_catch_flt1);
    wrt_r_complete_vector("Pred_ACindex", pred_Survindex);
    wrt_r_complete_vector("obs_ACindex", Acoustic_data);
    wrt_r_complete_vector("AC_yrs", Acoustic_yrs);
    wrt_r_complete_vector("Acoustic_sd", Acoustic_sd);
    wrt_r_complete_vector("SS_flt1", SS_flt1);
    wrt_r_complete_vector("SS_flt2", SS_flt2);
    wrt_r_complete_vector("b_ramp", b_recbias);
    wrt_r_complete_vector("pred_SurvBio", pred_SurvBio);
    wrt_r_complete_vector("pred_SurvBio_ind", pred_SurvBio_ind);
    open_r_matrix("Obs_fcomp");
      wrt_r_matrix(obs_fcomp);
    close_r_matrix();
    open_r_matrix("Pred_fcomp");
      wrt_r_matrix(pred_fcomp_compressed);
    close_r_matrix();
    open_r_matrix("Obs_Survcomp");
      wrt_r_matrix(obs_survcomp);
    close_r_matrix();
    open_r_matrix("Pred_Survcomp");
      wrt_r_matrix(pred_Survcomp_compressed);
    close_r_matrix();
    open_r_matrix("Fishery_Sel");
      wrt_r_matrix(Sel_fishery);
    close_r_matrix();
    wrt_r_complete_vector("Survey_Sel", Sel_surv);
    open_r_matrix("AE_means");
      wrt_r_matrix(AE_means);
    close_r_matrix();
    open_r_matrix("AE_SEs");
      wrt_r_matrix(AE_SEs);
    close_r_matrix();
    open_r_matrix("ObsComp_flt1");
      wrt_r_matrix(ObsComp_flt1);
    close_r_matrix();
    open_r_matrix("ObsComp_flt2");
      wrt_r_matrix(ObsComp_flt2);
    close_r_matrix();
    wrt_r_complete_vector("log_rec_devs", log_rec_devs);
    wrt_r_complete_vector("log_rec_devs_wbiasadj", log(exp_recdevs));
    wrt_r_complete_vector("log_Fint", log_Fint);
    wrt_r_complete_vector("fsel_par", fsel_par);
    wrt_r_complete_vector("survsel_par", survsel_par);
    open_r_matrix("sel_devs");
      wrt_r_matrix(sel_devs_mat);
    close_r_matrix();
    open_r_vector("params");
      wrt_r_item("log_survq",log_survq);          //Fishery Catchability 
      wrt_r_item("steepness",steepness);      //Steepness for recruitment
      wrt_r_item("log_R0",log_R0);          //Unfished recruitment for Florida cells
      wrt_r_item("log_ftheta",log_ftheta);           //log cv for fishery catch
      wrt_r_item("log_survtheta",log_survtheta);           //log cv for fishery catch
      wrt_r_item("log_surv_sd",log_surv_sd);            //log sd for recruitment
     close_r_vector();
    open_r_vector("NLL_Components");            //NLL components
      wrt_r_item("NLL_Catch",NLLcatch);     
      wrt_r_item("NLL_AgeComp",NLL1);      
      wrt_r_item("NLL_SurvIndex",NLL2);          
      wrt_r_item("NLL_SurvComp",NLL3);           
      wrt_r_item("NLL_M",NLL4);           
      wrt_r_item("NLL_Steep",NLL5);            
      wrt_r_item("NLL_Rec",NLL6);            
      wrt_r_item("NLL_SelDevs",NLL7);            
      wrt_r_item("NLL_DMFishery",NLL8);            
      wrt_r_item("NLL_DMSurv",NLL9);            
     close_r_vector();
  close_r_file();






