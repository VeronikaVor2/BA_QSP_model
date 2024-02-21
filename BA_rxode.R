## Description ----------------------------
## Physiological model of Bile acids distribution
## Author: Veronika Voronova
## published in [ Cell Mol Gastroenterol Hepatol. 2020;10(1):149-170]

## upload libraries ##
library(rxode2)
library(tidyverse)
theme_set(theme_bw(base_size = 16))
## compile the model ##
BA_model <- rxode2({
  # Model parameters ------------------
  
  ## Organ volumes
  PL   = 2.5  ; # [L] - plasma
  SIN  = 0.2  ; # [L] - sinusoidal space
  PV   = 0.45 ; # [L] - portal vein
  LIV  = 0.9  ; # [L] - liver
  BD   = 0.045; # [L] - bile duct
  GB   = 0.03 ; # [L] - gallbladder
  UINT = 0.2  ; # [L] - upper intestine
  LINT = 0.1  ; # [L] - lower intestine
  COL  = 0.3  ; # [L] - colon
  
  ## flows
  Q_MA       = 36   ; # [L/hour] - mesenteric artery plasma flow
  Q_PV       = 36   ; # [L/hour] - portal vein plasma flow
  Q_HA       = 12   ; # [L/hour] - hepatic artery plasma flow
  Q_HV       = 48   ; # [L/hour] - hepatic vein plasma flow
  Q_ex       = 10.8 ; # [L/hour] - biliary secretion (liver -> bile duct)
  Q_sec      = 0.027; # [L/hour] - biliary secretion (bile duct -> gall bladder)
  fr_fill_gb = 0.3  ; # [-] - BA fraction entering gall bladder
  
  ## transport constants
  kref    =   0  ; # [1/hour] - BA transport from liver to sinusoidal space
  kgbint  = 0.6  ; # [1/hour] - gall bladder emptying, depends on food type
  kintint = 0.18 ; # [1/hour] -  intestinal transit  
  kintcol = 0.12 ; # [1/hour] - intestinal transit  
  kcolel  = 0.07 ; # [1/hour] -  colonic transit 
  
  ## food-related changes
  gm      = 3   ; # [-] - gasric motolity increase in response to meal intake
  pc      = 0.5 ; # [-] - portal circulation increase in repsonse to meal intake
  tmbt    = 1   ; # [-] - microbiota activity modulation by food

  ftime_1 = 0.5 ; # [hours] - 1st food time
  ftime_2 = 4   ; # [hours] - 2nd food time
  ftime_3 = 10  ; # [hours] - 3rd food time
  ftime_4 = 12.5; #[hours]  - 4th food time
  
  ## BA-specific parameters 
  ksyn_ca     = 27.53  ; # [umol/hour] - CA synthesis rate
  ksyn_cdca   = 18.574 ; # [umol/hour] - CDCA synthesis rate
  kconj       = 4.8    ; # [1/hour] - BA conjugation
  fr_conj_tau = 0.25   ; # [-] - BA fraction conjugated with taurine
  
  klintdcj_tba = 1     ; # [1/hour] - taurine-conjugated BA deconjugation in the lower intestine
  klintdcj_gba = 4     ; # [1/hour] - glycine-conjugated BA deconjugation in the lower intestine
  kcoldcj      = 5     ; # [1/hour] - conjugated BA deconjugation in the colon
  
  khupt_uca     = 724  ; # [1/hour] - uCA hepatic uptake
  khupt_ucdca   = 243  ; # [1/hour] - uCDCA hepatic uptake
  khupt_udca    = 190  ; # [1/hour] - uDCA hepatic uptake
  khupt_cca     = 2016 ; # [1/hour] - tCA and gCA hepatic uptake
  khupt_ccdca   = 848  ; # [1/hour] - tCDCA and gCDCA hepatic uptake
  khupt_cdca    = 720  ; # [1/hour] - tDCA and gDCA hepatic uptake
  
  klintabs_c_cdcadca = 2.4                  ; # [1/hour] - tDCA, tCDCA, gDCA, gDCA absorption from the lower intestine
  klintabs_u_cdcadca = klintabs_c_cdcadca*3 ; # [1/hour] - uDCA and uCDCA absorption in the lower intestine
  klintabs_c_ca = klintabs_c_cdcadca*2.5    ; # [1/hour] - tCA and gCA absorption in the lower intestine
  klintabs_u_ca = klintabs_c_cdcadca*2.5    ; # [1/hour] - uCA absorption in the lower intestine
  
  kcolabs_cdca=klintabs_c_cdcadca*0.2       ; # [1/hour] - CDCA absorption from the colon
  kcolabs_ca=klintabs_c_cdcadca*0.02        ; # [1/hour] - CA absorption from the colon
  kcolabs_dca=klintabs_c_cdcadca*0.13       ; # [1/hour] - DCA absorption from the colon
  
  ## estimated
  kintabs_gba = 0.09     ; # [1/hour] - conjugated BA absorption from the upper intestine
  ksyn_dca    = 0.33     ; # [1/hour] - DCA synthesis rate (dehydroxylation)
  ksol_dca    = 0.041    ; # [1/hour] -  DCA solubilization rate 
  kel_cdca    = 0.6285   ; # [1/hour] - sCDCA loss due to secondary BA synthesis 
  lint_dcj    = 0.1373   ; # [1/hour] - BA deconjugation in lower intestine
  
  ## BA regulation parameters
  bln_c4         = 0.2   ; # [-] - auxilliary parameter to ensure daily AUC(C4)=1
  Fmax_fxr         = 2   ; # [-] - maximal FXR activation
  EC50_fxr          = 55   ; # [umol] - BA amount 
  kdel            = 0.4  ; # [1/hour] - delay between FXR activation and FGF-19 synthesis
  f_FXR_c_cdcadca = 1.7  ; # [-] - scaling factor for FXR activation by tDCA, tCDCA, gDCA, gDCA
  f_FXR_c_ca      = 5    ; # [-] - scaling factor for FXR activation by tCA, gCA
  f_FXR_u_cdca    = 4    ; # [-] - scaling factor for FXR activation by uCDCA
  f_FXR_u_dca     = 8    ; # [-] - scaling factor for FXR activation by uDCA
  a_ba_fgf19      = 2.8  ; # [-] - FGF-19 effect on BA synthesis
  ## equations --------------------------------------------------
  
  ## food intake impact 1 in postprandial state and 0 in fasting
  a1 = (0.5 < (ceil((t-ftime_1)/24)-ceil((t-ftime_1-3.5)/24)) + (ceil((t-ftime_2)/24)-ceil((t-ftime_2-3.5)/24))+ (ceil((t-ftime_3)/24)-ceil((t-ftime_3-3.5)/24))+ (ceil((t-ftime_4)/24)-ceil((t-ftime_4-3.5)/24)))
  a2 = (0.5 < (ceil((t-ftime_1)/24)-ceil((t-ftime_1-2)/24))+ (ceil((t-ftime_2)/24)-ceil((t-ftime_2-2)/24))+ (ceil((t-ftime_3)/24)-ceil((t-ftime_3-2)/24))+ (ceil((t-ftime_4)/24)-ceil((t-ftime_4-2)/24)))
  
  ## FXR activation 
  C4              = exp((1-fgf19)*a_ba_fgf19)  ## FCF-19 effect on C4 
  FXR_a_c_cdcadca = (tCDCA_lint+gCDCA_lint)+(tDCA_lint+gDCA_lint)/f_FXR_c_cdcadca
  FXR_a_c_ca      = (tCA_lint+gCA_lint)/f_FXR_c_ca 
  FXR_a_u_cdca    = CDCA_lint/f_FXR_u_cdca
  FXR_a_u_dca     = DCA_lint/f_FXR_u_dca 
  fxr_sum         = FXR_a_c_cdcadca + FXR_a_c_ca + FXR_a_u_cdca + FXR_a_u_dca
  FXRa            = Fmax_fxr*fxr_sum/EC50_fxr/(fxr_sum/EC50_fxr+1)
  
  ## control parameters
  gb_bsv       = 1  ; # [-] to simulate variability in gallbladder emptying rate
  intm_bsv     = 1  ; # [-] to simulate variability in intestinal motility
  colm_bsv     = 1  ; # [-] to simulate variability in colonic motility
  uint_abs_bsv = 1  ; # [-] to simulate variability in upper intestinal absorption
  lint_abs_bsv = 1  ; # [-] to simulate variability in lower intestinal absorption
  col_abs_bsv  = 1  ; # [-] to simulate variability in colonic absorption
  lint_dcj_bsv = 1  ; # [-] to simulate variability in lower intestinal deconjugation
  col_dcj_bsv  = 1  ; # [-] to simulate variability in colonic deconjugation
  mbt_bsv      = 1  ; # [-] to simulate variability in secondary BA production by microbiota
  hupt_bsv     = 1  ; # [-] to simulate variability in BA hepatic uptake
  ## concentrations calculation
  
  # individual BA
  CA_liv_um    = CA_liv/LIV	
  tCA_liv_um   = tCA_liv/LIV	
  gCA_liv_um   = gCA_liv/LIV
  CDCA_liv_um  = CDCA_liv/LIV	
  tCDCA_liv_um = tCDCA_liv/LIV	
  gCDCA_liv_um = gCDCA_liv/LIV	
  DCA_liv_um   = DCA_liv/LIV	
  tDCA_liv_um  = tDCA_liv/LIV	
  gDCA_liv_um  = gDCA_liv/LIV
  
  tCA_bd_um    = tCA_bd/BD	
  gCA_bd_um    = gCA_bd/BD	
  tCDCA_bd_um  = tCDCA_bd/BD	
  gCDCA_bd_um  = gCDCA_bd/BD	
  tDCA_bd_um   = tDCA_bd/BD	
  gDCA_bd_um   = gDCA_bd/BD	
  
  tCA_gb_um    = tCA_gb/GB	
  gCA_gb_um    = gCA_gb/GB	
  tCDCA_gb_um  = tCDCA_gb/GB	
  gCDCA_gb_um  = gCDCA_gb/GB	
  tDCA_gb_um   = tDCA_gb/GB	
  gDCA_gb_um   = gDCA_gb/GB
  
  tCA_uint_um   = tCA_uint/UINT	
  gCA_uint_um   = gCA_uint/UINT	
  tCDCA_uint_um = tCDCA_uint/UINT	
  gCDCA_uint_um = gCDCA_uint/UINT	
  tDCA_uint_um  = tDCA_uint/UINT	
  gDCA_uint_um  = gDCA_uint/UINT
  
  CA_lint_um    = CA_lint/LINT
  tCA_lint_um   = tCA_lint/LINT	
  gCA_lint_um   = gCA_lint/LINT
  CDCA_lint_um  = CDCA_lint/LINT	
  tCDCA_lint_um = tCDCA_lint/LINT	
  gCDCA_lint_um = gCDCA_lint/LINT	
  DCA_lint_um   = DCA_lint/LINT	
  tDCA_lint_um  = tDCA_lint/LINT	
  gDCA_lint_um  = gDCA_lint/LINT
  
  CA_pl_um    = CA_pl/PL	
  tCA_pl_um   = tCA_pl/PL	
  gCA_pl_um   = gCA_pl/PL
  CDCA_pl_um  = CDCA_pl/PL	
  tCDCA_pl_um = tCDCA_pl/PL	
  gCDCA_pl_um = gCDCA_pl/PL	
  DCA_pl_um   = DCA_pl/PL	
  tDCA_pl_um  = tDCA_pl/PL	
  gDCA_pl_um  = gDCA_pl/PL	
  
  CA_pv_um    = CA_pv/PV	
  tCA_pv_um   = tCA_pv/PV	
  gCA_pv_um   = gCA_pv/PV
  CDCA_pv_um  = CDCA_pv/PV	
  tCDCA_pv_um = tCDCA_pv/PV	
  gCDCA_pv_um = gCDCA_pv/PV	
  DCA_pv_um   = DCA_pv/PV	
  tDCA_pv_um  = tDCA_pv/PV	
  gDCA_pv_um  = gDCA_pv/PV	
  
  CA_col_um     = CA_col/COL	
  tCA_col_um    = tCA_col/COL
  gCA_col_um    = gCA_col/COL
  CDCA_col_um   = CDCA_col/COL	
  tCDCA_col_um  = tCDCA_col/COL	
  gCDCA_col_um  = gCDCA_col/COL	
  DCA_col_um    = DCA_col/COL	
  tDCA_col_um   = tDCA_col/COL	
  gDCA_col_um   = gDCA_col/COL
  DCAins_col_um = DCAins_col/COL
  
  # total BA
  totCA_pl_um     = CA_pl_um+tCA_pl_um+gCA_pl_um
  totCDCA_pl_um   = CDCA_pl_um+tCDCA_pl_um+gCDCA_pl_um
  totDCA_pl_um    = DCA_pl_um+tDCA_pl_um+gDCA_pl_um
  totCA_pv_um     = CA_pv_um+tCA_pv_um+gCA_pv_um
  totCDCA_pv_um   = CDCA_pv_um+tCDCA_pv_um+gCDCA_pv_um
  totDCA_pv_um    = DCA_pv_um+tDCA_pv_um+gDCA_pv_um
  totCA_liv_um    = CA_liv_um+tCA_liv_um+gCA_liv_um
  totCDCA_liv_um  = CDCA_liv_um+tCDCA_liv_um+gCDCA_liv_um
  totDCA_liv_um   = DCA_liv_um+tDCA_liv_um+gDCA_liv_um
  totCA_gb_um     = tCA_gb_um+gCA_gb_um
  totCDCA_gb_um   = tCDCA_gb_um+gCDCA_gb_um
  totDCA_gb_um    = tDCA_gb_um+gDCA_gb_um
  totCA_uint_um   = tCA_uint_um+gCA_uint_um
  totCDCA_uint_um = tCDCA_uint_um+gCDCA_uint_um
  totDCA_uint_um  = tDCA_uint_um+gDCA_uint_um
  totCA_lint_um   = CA_lint_um+tCA_lint_um+gCA_lint_um
  totCDCA_lint_um = CDCA_lint_um+tCDCA_lint_um+gCDCA_lint_um
  totDCA_lint_um  = DCA_lint_um+tDCA_lint_um+gDCA_lint_um
  totCA_col_um    = CA_col_um+tCA_col_um+gCA_col_um
  totCDCA_col_um  = CDCA_col_um+tCDCA_col_um+gCDCA_col_um
  totDCA_col_um   = DCA_col_um+tDCA_col_um+gDCA_col_um ; +DCAins_col_um
  
  totTAU_pl    = tCA_pl_um+tCDCA_pl_um+tDCA_pl_um
  totGLY_pl    = gCA_pl_um+gCDCA_pl_um+gDCA_pl_um
  totUN_pl     = CA_pl_um+CDCA_pl_um+DCA_pl_um
  totTAU_pv    = tCA_pv_um+tCDCA_pv_um+tDCA_pv_um
  totGLY_pv    = gCA_pv_um+gCDCA_pv_um+gDCA_pv_um
  totUN_pv     = CA_pv_um+CDCA_pv_um+DCA_pv_um
  totTAU_col   = tCA_col_um+tCDCA_col_um+tDCA_col_um
  totGLY_col   = gCA_col_um+gCDCA_col_um+gDCA_col_um
  totUN_col    = CA_col_um+CDCA_col_um+DCA_col_um+DCAins_col_um
  totTAU_lint  = tCA_lint_um+tCDCA_lint_um+tDCA_lint_um
  totGLY_lint  = gCA_lint_um+gCDCA_lint_um+gDCA_lint_um
  totUN_lint   = CA_lint_um+CDCA_lint_um+DCA_lint_um
  
  TBA_pl_um    = CA_pl_um+tCA_pl_um+gCA_pl_um+CDCA_pl_um+tCDCA_pl_um+gCDCA_pl_um+DCA_pl_um+tDCA_pl_um+gDCA_pl_um
  TBA_pv_um    = CA_pv_um+tCA_pv_um+gCA_pv_um+CDCA_pv_um+tCDCA_pv_um+gCDCA_pv_um+DCA_pv_um+tDCA_pv_um+gDCA_pv_um
  TBA_liv_um   = CA_liv_um+tCA_liv_um+gCA_liv_um+CDCA_liv_um+tCDCA_liv_um+gCDCA_liv_um+DCA_liv_um+tDCA_liv_um+gDCA_liv_um
  TBA_gb_um    = tCA_gb_um+gCA_gb_um+tCDCA_gb_um+gCDCA_gb_um+tDCA_gb_um+gDCA_gb_um
  TBA_bd_um    = tCA_bd_um+gCA_bd_um+tCDCA_bd_um+gCDCA_bd_um+tDCA_bd_um+gDCA_bd_um
  TBA_col_um   = (CA_col_um+tCA_col_um+gCA_col_um+CDCA_col_um+tCDCA_col_um+gCDCA_col_um+DCA_col_um+tDCA_col_um+gDCA_col_um+DCAins_col_um)
  TBA_uint_um  = (tCA_uint+gCA_uint+tCDCA_uint+gCDCA_uint+tDCA_uint+gDCA_uint)/UINT
  TBA_lint_um  = (CA_lint+tCA_lint+gCA_lint+CDCA_lint+tCDCA_lint+gCDCA_lint+DCA_lint+tDCA_lint+gDCA_lint)/LINT

  TBA_pl_bln_ch = TBA_pl_um/PL
  
  ## BA percentage calculation 
  CA_pl_perc     = (CA_pl+tCA_pl+gCA_pl)/PL/(TBA_pl_um)*100
  CDCA_pl_perc   = (CDCA_pl+tCDCA_pl+gCDCA_pl)/PL/(TBA_pl_um)*100
  DCA_pl_perc    = (DCA_pl+tDCA_pl+gDCA_pl)/PL/(TBA_pl_um)*100
  CA_liv_perc    = (CA_liv+tCA_liv+gCA_liv)/LIV/(TBA_liv_um)*100
  CDCA_liv_perc  = (CDCA_liv+tCDCA_liv+gCDCA_liv)/LIV/(TBA_liv_um)*100
  DCA_liv_perc   = (DCA_liv+tDCA_liv+gDCA_liv)/LIV/(TBA_liv_um)*100
  CA_pv_perc     = (CA_pv+tCA_pv+gCA_pv)/PV/(TBA_pv_um)*100
  CDCA_pv_perc   = (CDCA_pv+tCDCA_pv+gCDCA_pv)/PV/(TBA_pv_um)*100
  DCA_pv_perc    = (DCA_pv+tDCA_pv+gDCA_pv)/PV/(TBA_pv_um)*100
  CA_uint_perc   = (tCA_uint+gCA_uint)/UINT/(TBA_uint_um)*100
  CDCA_uint_perc = (tCDCA_uint+gCDCA_uint)/UINT/(TBA_uint_um)*100
  DCA_uint_perc  = (tDCA_uint+gDCA_uint)/UINT/(TBA_uint_um)*100
  CA_lint_perc   = (CA_lint+tCA_lint+gCA_lint)/LINT/(TBA_lint_um)*100
  CDCA_lint_perc = (CDCA_lint+tCDCA_lint+gCDCA_lint)/LINT/(TBA_lint_um)*100
  DCA_lint_perc  = (DCA_lint+tDCA_lint+gDCA_lint)/LINT/(TBA_lint_um)*100	
  CA_col_perc    = (CA_col+tCA_col+gCA_col)/COL/(TBA_col_um)*100
  CDCA_col_perc  = (CDCA_col+tCDCA_col+gCDCA_col)/COL/(TBA_col_um)*100
  DCA_col_perc   = (DCA_col+tDCA_col+gDCA_col+DCAins_col)/COL/(TBA_col_um)*100	
  CA_gb_perc     = (tCA_gb+gCA_gb)/GB/(TBA_gb_um)*100
  CDCA_gb_perc   = (tCDCA_gb+gCDCA_gb)/GB/(TBA_gb_um)*100
  DCA_gb_perc    = (tDCA_gb+gDCA_gb)/GB/(TBA_gb_um)*100
  
  tauBA_pl_perc   = (tCA_pl_um+tCDCA_pl_um+tDCA_pl_um)/TBA_pl_um*100
  gBA_pl_perc     = (gCA_pl_um+gCDCA_pl_um+gDCA_pl_um)/TBA_pl_um*100
  uBA_pl_perc     = (CA_pl_um+CDCA_pl_um+DCA_pl_um)/TBA_pl_um*100
  tauBA_pv_perc   = (tCA_pv_um+tCDCA_pv_um+tDCA_pv_um)/TBA_pv_um*100
  gBA_pv_perc     = (gCA_pv_um+gCDCA_pv_um+gDCA_pv_um)/TBA_pv_um*100
  uBA_pv_perc     = (CA_pv_um+CDCA_pv_um+DCA_pv_um)/TBA_pv_um*100
  tauBA_lint_perc = (tCA_lint+tCDCA_lint+tDCA_lint)/LINT/TBA_lint_um*100
  gBA_lint_perc   = (gCA_lint+gCDCA_lint+gDCA_lint)/LINT/TBA_lint_um*100
  uBA_lint_perc   = (CA_lint+CDCA_lint+DCA_lint)/LINT/TBA_lint_um*100
  tauBA_col_perc  = (tCA_col_um+tCDCA_col_um+tDCA_col_um)/TBA_col_um*100
  gBA_col_perc    = (gCA_col_um+gCDCA_col_um+gDCA_col_um)/TBA_col_um*100
  uBA_col_perc    = (CA_col_um+CDCA_col_um+DCA_col_um+DCAins_col_um)/TBA_col_um*100
  tauBA_gb_perc   = (tCA_gb+tCDCA_gb+tDCA_gb)/GB/(TBA_gb_um)*100
  gBA_gb_perc     = (gCA_gb+gCDCA_gb+gDCA_gb)/GB/(TBA_gb_um)*100
  tauBA_uint_perc = (tCA_uint+tCDCA_uint+tDCA_uint)/UINT/(TBA_uint_um)*100
  gBA_uint_perc   = (gCA_uint+gCDCA_uint+gDCA_uint)/UINT/(TBA_uint_um)*100
  
  cBA_pl_perc     = tauBA_pl_perc+gBA_pl_perc
  cBA_pv_perc     = tauBA_pv_perc+gBA_pv_perc
  cBA_uint_perc   = tauBA_uint_perc+gBA_uint_perc
  cBA_lint_perc   = tauBA_lint_perc+gBA_lint_perc
  cBA_col_perc    = tauBA_col_perc+gBA_col_perc
  cBA_gb_perc     = tauBA_gb_perc+gBA_gb_perc
  
  ## Reaction rates ------------------------------
  
  # synthesis and conjugation 
  CAsyn        = ksyn_ca*(C4 + bln_c4)
  CDCAsyn      = ksyn_cdca
  
  tCAconj      = (kconj*fr_conj_tau)*CA_liv
  tCDCAconj    = (kconj*fr_conj_tau)*CDCA_liv
  tDCAconj     = (kconj*fr_conj_tau)*DCA_liv
  
  gCAconj      = (kconj*(1-fr_conj_tau))*CA_liv
  gCDCAconj    = (kconj*(1-fr_conj_tau))*CDCA_liv
  gDCAconj     = (kconj*(1-fr_conj_tau))*DCA_liv
  
  # secretion 
  tCAseq_tot   = (Q_ex/LIV)*tCA_liv
  tCDCAseq_tot = (Q_ex/LIV)*tCDCA_liv
  tDCAseq_tot  = (Q_ex/LIV)*tDCA_liv
  
  gCAseq_tot   = (Q_ex/LIV)*gCA_liv
  gCDCAseq_tot = (Q_ex/LIV)*gCDCA_liv
  gDCAseq_tot  = (Q_ex/LIV)*gDCA_liv
  
  tCAgb_fill   = ((Q_sec/BD)*fr_fill_gb)*tCA_bd
  tCDCAgb_fill = ((Q_sec/BD)*fr_fill_gb)*tCDCA_bd
  tDCAgb_fill  = ((Q_sec/BD)*fr_fill_gb)*tDCA_bd
  
  gCAgb_fill   = ((Q_sec/BD)*fr_fill_gb)*gCA_bd
  gCDCAgb_fill = ((Q_sec/BD)*fr_fill_gb)*gCDCA_bd
  gDCAgb_fill  = ((Q_sec/BD)*fr_fill_gb)*gDCA_bd
  
  tCAseqf      = ((Q_sec/BD)*(1-fr_fill_gb))*tCA_bd
  tCDCAseqf    = ((Q_sec/BD)*(1-fr_fill_gb))*tCDCA_bd
  tDCAseqf     = ((Q_sec/BD)*(1-fr_fill_gb))*tDCA_bd
  
  gCAseqf      = ((Q_sec/BD)*(1-fr_fill_gb))*gCA_bd
  gCDCAseqf    = ((Q_sec/BD)*(1-fr_fill_gb))*gCDCA_bd
  gDCAseqf     = ((Q_sec/BD)*(1-fr_fill_gb))*gDCA_bd
  
  tCAseqpp     = kgbint*tCA_gb*a2*gb_bsv
  tCDCAseqpp   = kgbint*tCDCA_gb*a2*gb_bsv
  tDCAseqpp    = kgbint*tDCA_gb*a2*gb_bsv
  
  gCAseqpp     = kgbint*gCA_gb*a2*gb_bsv
  gCDCAseqpp   = kgbint*gCDCA_gb*a2*gb_bsv
  gDCAseqpp    = kgbint*gDCA_gb*a2*gb_bsv
  
  # gastro-intestinal fluxes
  tCAintint    = kintint*tCA_uint*(1+a1*gm)*intm_bsv
  tCDCAintint  = kintint*tCDCA_uint*(1+a1*gm)*intm_bsv
  tDCAintint   = kintint*tDCA_uint*(1+a1*gm)*intm_bsv
  
  gCAintint    = kintint*gCA_uint*(1+a1*gm)*intm_bsv
  gCDCAintint  = kintint*gCDCA_uint*(1+a1*gm)*intm_bsv
  gDCAintint   = kintint*gDCA_uint*(1+a1*gm)*intm_bsv
  
  tCAintcol    = kintcol*tCA_lint*(1+a1*gm)*intm_bsv
  tCDCAintcol  = kintcol*tCDCA_lint*(1+a1*gm)*intm_bsv
  tDCAintcol   = kintcol*tDCA_lint*(1+a1*gm)*intm_bsv
  
  gCAintcol    = kintcol*gCA_lint*(1+a1*gm)*intm_bsv
  gCDCAintcol  = kintcol*gCDCA_lint*(1+a1*gm)*intm_bsv
  gDCAintcol   = kintcol*gDCA_lint*(1+a1*gm)*intm_bsv
  
  CAintcol     = kintcol*CA_lint*(1+a1*gm)*intm_bsv
  CDCAintcol   = kintcol*CDCA_lint*(1+a1*gm)*intm_bsv
  DCAintcol    = kintcol*DCA_lint*(1+a1*gm)*intm_bsv
  
  CAfec        = kcolel*CA_col*colm_bsv
  CDCAfec      = kcolel*CDCA_col*colm_bsv
  DCAfec       = kcolel*DCA_col*colm_bsv
  DCAinsfec    = kcolel*DCAins_col*colm_bsv
  
  # reabsorption
  gCDCAintabs  = kintabs_gba*gCDCA_uint*uint_abs_bsv
  gDCAintabs   =  kintabs_gba*gDCA_uint*uint_abs_bsv
  
  gCAlintabs   = klintabs_c_ca*gCA_lint *lint_abs_bsv
  gCDCAlintabs = klintabs_c_cdcadca*gCDCA_lint *lint_abs_bsv
  gDCAlintabs  =  klintabs_c_cdcadca*gDCA_lint*lint_abs_bsv
  
  tCAlintabs   = klintabs_c_ca*tCA_lint*lint_abs_bsv
  tCDCAlintabs = klintabs_c_cdcadca*tCDCA_lint*lint_abs_bsv
  tDCAlintabs  = klintabs_c_cdcadca*tDCA_lint*lint_abs_bsv
  
  CAlintabs    = klintabs_u_ca*CA_lint*lint_abs_bsv
  CDCAlintabs  = klintabs_u_cdcadca*CDCA_lint*lint_abs_bsv
  DCAlintabs   =  klintabs_u_cdcadca*DCA_lint*lint_abs_bsv
  
  CAcolabs     =  kcolabs_ca*CA_col*col_abs_bsv
  CDCAcolabs   =  kcolabs_cdca*CDCA_col*col_abs_bsv
  DCAcolabs    =  kcolabs_dca*DCA_col*col_abs_bsv
  
  # deconjugation and dehydroxylation
  tCAdcjlint   = klintdcj_tba*tCA_lint *lint_dcj*(1-a1*tmbt)*lint_dcj_bsv
  tCDCAdcjlint = klintdcj_tba*tCDCA_lint *lint_dcj*(1-a1*tmbt)*lint_dcj_bsv 
  tDCAdcjlint  = klintdcj_tba*tDCA_lint *lint_dcj*(1-a1*tmbt)*lint_dcj_bsv
  
  gCAdcjlint   = klintdcj_gba*gCA_lint*lint_dcj*(1-a1*tmbt)*lint_dcj_bsv
  gCDCAdcjlint = klintdcj_gba*gCDCA_lint *lint_dcj*(1-a1*tmbt)*lint_dcj_bsv
  gDCAdcjlint  = klintdcj_gba*gDCA_lint *lint_dcj*(1-a1*tmbt)*lint_dcj_bsv
  
  tCAdcjcol    = kcoldcj*tCA_col*col_dcj_bsv
  tCDCAdcjcol  = kcoldcj*tCDCA_col*col_dcj_bsv
  tDCAdcjcol   = kcoldcj*tDCA_col*col_dcj_bsv
  
  gCAdcjcol    = kcoldcj*gCA_col*col_dcj_bsv
  gCDCAdcjcol  = kcoldcj*gCDCA_col*col_dcj_bsv
  gDCAdcjcol   = kcoldcj*gDCA_col*col_dcj_bsv
  
  DCAsyn       = ksyn_dca*CA_col*mbt_bsv
  DCAsol       = ksol_dca*DCAins_col
  CDCAelim     = kel_cdca*CDCA_col*mbt_bsv
  
  # plasma distribution
  tCApvsin      = (Q_PV/PV)*tCA_pv*(1+a1*pc)
  tCDCApvsin   = (Q_PV/PV)*tCDCA_pv*(1+a1*pc)
  tDCApvsin    = (Q_PV/PV)*tDCA_pv*(1+a1*pc)
  
  gCApvsin     = (Q_PV/PV)*gCA_pv*(1+a1*pc)
  gCDCApvsin   = (Q_PV/PV)*gCDCA_pv*(1+a1*pc)
  gDCApvsin    = (Q_PV/PV)*gDCA_pv*(1+a1*pc)
  
  CApvsin      = (Q_PV/PV)*CA_pv*(1+a1*pc)
  CDCApvsin    = (Q_PV/PV)*CDCA_pv*(1+a1*pc)
  DCApvsin     = (Q_PV/PV)*DCA_pv*(1+a1*pc)
  
  tCAsinpl     = (Q_HV/SIN)*tCA_sin*(1+a1*pc) - (Q_HA/PL)*tCA_pl*(1+a1*pc)
  tCDCAsinpl   = (Q_HV/SIN)*tCDCA_sin*(1+a1*pc) - (Q_HA/PL)*tCDCA_pl*(1+a1*pc)
  tDCAsinpl    = (Q_HV/SIN)*tDCA_sin*(1+a1*pc) - (Q_HA/PL)*tDCA_pl*(1+a1*pc)
  
  gCAsinpl     = (Q_HV/SIN)*gCA_sin*(1+a1*pc) - (Q_HA/PL)*gCA_pl*(1+a1*pc)
  gCDCAsinpl   = (Q_HV/SIN)*gCDCA_sin*(1+a1*pc) - (Q_HA/PL)*gCDCA_pl*(1+a1*pc)
  gDCAsinpl    = (Q_HV/SIN)*gDCA_sin*(1+a1*pc) - (Q_HA/PL)*gDCA_pl*(1+a1*pc)
  
  CAsinpl      = (Q_HV/SIN)*CA_sin*(1+a1*pc) - (Q_HA/PL)*CA_pl*(1+a1*pc)
  CDCAsinpl    = (Q_HV/SIN)*CDCA_sin*(1+a1*pc) - (Q_HA/PL)*CDCA_pl*(1+a1*pc)
  DCAsinpl     = (Q_HV/SIN)*DCA_sin*(1+a1*pc) - (Q_HA/PL)*DCA_pl*(1+a1*pc)
  
  tCAplpv      = (Q_MA/PL)*tCA_pl*(1+a1*pc)
  tCDCAplpv    = (Q_MA/PL)*tCDCA_pl*(1+a1*pc)
  tDCAplpv     = (Q_MA/PL)*tDCA_pl*(1+a1*pc)
  
  # hepatic uptake
  tCAhupt      = khupt_cca*tCA_sin*hupt_bsv - kref*tCA_liv
  gCAhupt      = khupt_cca*gCA_sin*hupt_bsv - kref*gCA_liv
  CAhupt       = khupt_uca*CA_sin - kref*CA_liv
  
  tCDCAhupt    = khupt_ccdca*tCDCA_sin*hupt_bsv - kref*tCDCA_liv
  gCDCAhupt    = khupt_ccdca*gCDCA_sin*hupt_bsv - kref*gCDCA_liv
  CDCAhupt     = khupt_ucdca*CDCA_sin - kref*CDCA_liv
  
  tDCAhupt     = khupt_cdca*tDCA_sin*hupt_bsv - kref*tDCA_liv
  gDCAhupt     = khupt_cdca*gDCA_sin*hupt_bsv - kref*gDCA_liv
  DCAhupt      = khupt_udca*DCA_sin - kref*DCA_liv
  gCAplpv      = (Q_MA/PL)*gCA_pl*(1+a1*pc)
  gCDCAplpv    = (Q_MA/PL)*gCDCA_pl*(1+a1*pc)
  gDCAplpv     = (Q_MA/PL)*gDCA_pl*(1+a1*pc)
  
  CAplpv       = (Q_MA/PL)*CA_pl*(1+a1*pc)
  CDCAplpv     = (Q_MA/PL)*CDCA_pl*(1+a1*pc)
  DCAplpv      = (Q_MA/PL)*DCA_pl*(1+a1*pc)
  
  ## ordinary differential equation system --------------------------------
  # liver
  d/dt(CA_liv)    =   CAsyn +CAhupt - tCAconj - gCAconj
  d/dt(CDCA_liv)  = CDCAsyn + CDCAhupt  - tCDCAconj - gCDCAconj
  d/dt(DCA_liv)   = DCAhupt  - tDCAconj - gDCAconj
  
  d/dt(tCA_liv)   = tCAconj + tCAhupt- tCAseq_tot
  d/dt(tCDCA_liv) = tCDCAconj + tCDCAhupt- tCDCAseq_tot
  d/dt(tDCA_liv)  = tDCAconj + tDCAhupt- tDCAseq_tot
  
  d/dt(gCA_liv)   = gCAconj + gCAhupt- gCAseq_tot
  d/dt(gCDCA_liv) = gCDCAconj + gCDCAhupt- gCDCAseq_tot
  d/dt(gDCA_liv)  = gDCAconj + gDCAhupt- gDCAseq_tot
  
  # billiary tract
  d/dt(tCA_bd)    = tCAseq_tot - tCAgb_fill - tCAseqf
  d/dt(tCDCA_bd)  = tCDCAseq_tot - tCDCAgb_fill - tCDCAseqf
  d/dt(tDCA_bd)   = tDCAseq_tot - tDCAgb_fill - tDCAseqf
  
  d/dt(gCA_bd)    = gCAseq_tot - gCAgb_fill - gCAseqf
  d/dt(gCDCA_bd)  = gCDCAseq_tot - gCDCAgb_fill - gCDCAseqf
  d/dt(gDCA_bd)   = gDCAseq_tot - gDCAgb_fill - gDCAseqf
  
  d/dt(tCA_gb)    = tCAgb_fill - tCAseqpp
  d/dt(tCDCA_gb)  = tCDCAgb_fill - tCDCAseqpp
  d/dt(tDCA_gb)   = tDCAgb_fill - tDCAseqpp
  
  d/dt(gCA_gb)    = gCAgb_fill - gCAseqpp
  d/dt(gCDCA_gb)  = gCDCAgb_fill - gCDCAseqpp
  d/dt(gDCA_gb)   = gDCAgb_fill - gDCAseqpp
  
  # gastro-intestinal tract
  d/dt(tCA_uint)   = tCAseqpp + tCAseqf - tCAintint
  d/dt(tCDCA_uint) = tCDCAseqpp + tCDCAseqf - tCDCAintint
  d/dt(tDCA_uint)  = tDCAseqpp + tDCAseqf - tDCAintint
  
  d/dt(gCA_uint)   = gCAseqpp + gCAseqf - gCAintint
  d/dt(gCDCA_uint) = gCDCAseqpp + gCDCAseqf - gCDCAintint -gCDCAintabs
  d/dt(gDCA_uint)  = gDCAseqpp + gDCAseqf - gDCAintint -gDCAintabs
  
  d/dt(tCA_lint)   = tCAintint - tCAintcol - tCAlintabs - tCAdcjlint
  d/dt(tCDCA_lint) = tCDCAintint - tCDCAintcol - tCDCAlintabs - tCDCAdcjlint
  d/dt(tDCA_lint)  = tDCAintint - tDCAintcol - tDCAlintabs - tDCAdcjlint
  
  d/dt(gCA_lint)   = gCAintint - gCAintcol - gCAlintabs - gCAdcjlint
  d/dt(gCDCA_lint) = gCDCAintint - gCDCAintcol - gCDCAlintabs - gCDCAdcjlint
  d/dt(gDCA_lint)  = gDCAintint - gDCAintcol - gDCAlintabs - gDCAdcjlint
  
  d/dt(CA_lint)    = tCAdcjlint + gCAdcjlint - CAintcol - CAlintabs
  d/dt(CDCA_lint)  = tCDCAdcjlint + gCDCAdcjlint - CDCAintcol - CDCAlintabs
  d/dt(DCA_lint)   = tDCAdcjlint + gDCAdcjlint - DCAintcol - DCAlintabs
  
  d/dt(tCA_col)    = tCAintcol - tCAdcjcol 
  d/dt(tCDCA_col)  = tCDCAintcol - tCDCAdcjcol
  d/dt(tDCA_col)   = tDCAintcol - tDCAdcjcol
  
  d/dt(gCA_col)    = gCAintcol - gCAdcjcol 
  d/dt(gCDCA_col)  = gCDCAintcol - gCDCAdcjcol
  d/dt(gDCA_col)   = gDCAintcol - gDCAdcjcol
  
  d/dt(CA_col)     = CAintcol + tCAdcjcol + gCAdcjcol - DCAsyn- CAfec - CAcolabs
  d/dt(CDCA_col)   = CDCAintcol + tCDCAdcjcol + gCDCAdcjcol - CDCAelim - CDCAfec - CDCAcolabs
  d/dt(DCAins_col) = DCAsyn - DCAsol - DCAinsfec
  d/dt(DCA_col)    = DCAintcol + tDCAdcjcol + gDCAdcjcol + DCAsol - DCAfec - DCAcolabs
  
  # circulation
  d/dt(tCA_pv)     = tCAplpv+ tCAlintabs - tCApvsin
  d/dt(tCDCA_pv)   = tCDCAplpv+ tCDCAlintabs - tCDCApvsin
  d/dt(tDCA_pv)    = tDCAplpv+ tDCAlintabs - tDCApvsin
  
  d/dt(gCA_pv)     = gCAplpv+ gCAlintabs - gCApvsin
  d/dt(gCDCA_pv)   = gCDCAplpv+ gCDCAintabs+ gCDCAlintabs - gCDCApvsin
  d/dt(gDCA_pv)    = gDCAplpv+ gDCAintabs + gDCAlintabs - gDCApvsin
  
  d/dt(CA_pv)      = CAplpv + CAlintabs + CAcolabs - CApvsin
  d/dt(CDCA_pv)    = CDCAplpv + CDCAlintabs + CDCAcolabs - CDCApvsin
  d/dt(DCA_pv)     = DCAplpv + DCAlintabs + DCAcolabs - DCApvsin
  
  d/dt(tCA_sin)    = tCApvsin - tCAsinpl - tCAhupt
  d/dt(tCDCA_sin)  = tCDCApvsin - tCDCAsinpl - tCDCAhupt
  d/dt(tDCA_sin)   = tDCApvsin - tDCAsinpl - tDCAhupt
  
  d/dt(gCA_sin)    = gCApvsin - gCAsinpl - gCAhupt
  d/dt(gCDCA_sin)  = gCDCApvsin - gCDCAsinpl - gCDCAhupt
  d/dt(gDCA_sin)   = gDCApvsin - gDCAsinpl - gDCAhupt
  
  d/dt(CA_sin)     = CApvsin - CAsinpl - CAhupt
  d/dt(CDCA_sin)   = CDCApvsin - CDCAsinpl - CDCAhupt
  d/dt(DCA_sin)    = DCApvsin - DCAsinpl - DCAhupt
  
  d/dt(tCA_pl)     = tCAsinpl - tCAplpv
  d/dt(tCDCA_pl)   = tCDCAsinpl - tCDCAplpv
  d/dt(tDCA_pl)    = tDCAsinpl - tDCAplpv
  
  d/dt(gCA_pl)     = gCAsinpl - gCAplpv
  d/dt(gCDCA_pl)   = gCDCAsinpl - gCDCAplpv
  d/dt(gDCA_pl)    = gDCAsinpl - gDCAplpv
  
  d/dt(CA_pl)      = CAsinpl - CAplpv
  d/dt(CDCA_pl)    = CDCAsinpl - CDCAplpv
  d/dt(DCA_pl)     = DCAsinpl - DCAplpv
 
  # regulation
  d/dt(fgf19)      = kdel*(FXRa-fgf19)
})

### model-based simulations ---------------------------------------


## basic simulation
basic_evtab <- data.frame(Time = seq(0, 1680,0.01))
basic_sim   <- rxSolve(BA_model,basic_evtab)

## check food control functions
ggplot(basic_sim %>% filter(time > 1680-24)) + geom_line(aes(x = time, y = a1))
ggplot(basic_sim %>% filter(time > 1680-24)) + geom_line(aes(x = time, y = a2))

## check BA dynamics reproduction (fig. 3 from the manuscript)
vname_ba <- c('CA_liv_um','tCA_liv_um','gCA_liv_um','CDCA_liv_um','tCDCA_liv_um','gCDCA_liv_um','DCA_liv_um','tDCA_liv_um','gDCA_liv_um',
             'CA_pl_um','tCA_pl_um','gCA_pl_um','CDCA_pl_um','tCDCA_pl_um','gCDCA_pl_um','DCA_pl_um','tDCA_pl_um','gDCA_pl_um',
             'CA_pv_um','tCA_pv_um','gCA_pv_um','CDCA_pv_um','tCDCA_pv_um','gCDCA_pv_um','DCA_pv_um','tDCA_pv_um','gDCA_pv_um',
             'tCA_uint_um','gCA_uint_um','tCDCA_uint_um','gCDCA_uint_um','tDCA_uint_um','gDCA_uint_um',
             'CA_lint_um','tCA_lint_um','gCA_lint_um','CDCA_lint_um','tCDCA_lint_um','gCDCA_lint_um','DCA_lint_um','tDCA_lint_um','gDCA_lint_um',
             'CA_col_um','tCA_col_um','gCA_col_um','CDCA_col_um','tCDCA_col_um','gCDCA_col_um','DCA_col_um','tDCA_col_um','gDCA_col_um',
             'tCA_gb_um','gCA_gb_um','tCDCA_gb_um','gCDCA_gb_um','tDCA_gb_um','gDCA_gb_um','fgf19','C4')

basic_sim_ba <- basic_sim %>% 
  gather('variable','value',-time) %>% 
  filter(variable %in% vname_ba & time > max(time - 24)) %>% 
  separate(variable, c('variable','comp','unit'), sep = '_')

ggplot(basic_sim_ba %>% filter(!variable %in% c('fgf19','C4'))) +
  facet_wrap(~comp, scales = 'free') +
  geom_area(aes(x = time, y = value, fill = variable, group = variable), 
            color = 'black', alpha = 0.5) +
  xlab('Time, hours') + ylab('Concentration, uM')

##  check C4 and FGF19 dynamics
ggplot(basic_sim_ba %>% filter(variable %in% c('fgf19','C4'))) +
  geom_line(aes(x = time, y = value, color = variable, group = variable), size = 1) +
  xlab('Time, hours') + ylab('Relative change, -')


