*///////////////////////// STAGE 1 - data collection /////////////////////////;
data proj1;
   set '/home/u63638642/sasuser.v94/project1/data/z02.sas7bdat';
run;

*Investigating distribution of missing data;
data lb05;
 set proj1;
 subjid = _N_;
 if trt=1 then trt2="A"; *experimental;
 if trt=2 then trt2="P"; *control;
 chg=y_t2-y_t1;
 if y_t2="" then m=1;
 else m=0;
 if m=1 then chg2=-16;
 else chg2=chg;
run;

*descriptive;
proc means data=lb05;
run;

*///////////////////////// STAGE 2 - initial analysis of the raw data /////////////////////////;
*graphical analysis;

 proc sgplot data=lb05;
  styleattrs datacontrastcolors=(blue red) datasymbols=(circlefilled );
  scatter x=y_t1 y=chg2 / group=m;
 run;
 proc sgplot data=lb05;
  styleattrs datacontrastcolors=(blue red) datasymbols=(circlefilled );
  scatter x=y_t1 y=chg2/ group=trt2;
 run;
  proc sgplot data=lb05;
  styleattrs datacontrastcolors=(blue red) datasymbols=(circlefilled );
  scatter x=y_t1 y=chg2/ group=x2;
 run;
 
 
*///////////////////////// STAGE 3 - check for missing data distribution /////////////////////////; 
*we want to check distribution of missing data in each of the treatment categories;
proc sql;
    create table m_test as 
    select m, trt2 from lb05;
quit;

proc freq data=m_test; *here we have freq table for missings;
tables trt2*m / crosslist;
run;

proc means data=lb05;
run;


*///////////////////////// STAGE 4 - investigation of distribution of data /////////////////////////;
*we check what distribution of Y_T2 we have;
PROC UNIVARIATE DATA=lb05; 
var Y_T2; 
HISTOGRAM / normal; 
RUN;


*///////////////////////// STAGE 5 - logistic regression for tagret determining misingness /////////////////////////;
*logistic regression of missings;
proc logistic data=lb05;
  class trt2 x2;
  model m(Event='1') = trt2 y_t1 x1 x2 / link=logit scale=none;
  output out=LogisticOutput predicted=Probability;
run;


*///////////////////////// STAGE 6 - MCAR test /////////////////////////;
*MCAR;
proc mixed data=lb05;
 class trt2;
 model chg = y_t1 trt2 / solution;
 lsmeans trt2 / pdiff=all adjust=tukey;
run;


*///////////////////////// STAGE 7 - MAR test /////////////////////////;
*MAR - ORIGINAL;
proc mi data=lb05 out=out01 nimpute=500 seed=2012;
 class trt2 x2; 
 var trt2 y_t1 x1 x2 y_t2;
 monotone reg (y_t2 / details); *for innovation change here;
run;

*calculate chg;
data out02;
 set out01;
 chg = y_t2-y_t1;
run;

proc sql; *saved for later;
	create table fin_imput as
	select y_t1, x1, trt, x2, y_t2, trt2, chg, m
	from out02
	where _Imputation_=1;*randomly selected;
run;

*estimate model;
ods select none;
proc mixed data=out02;
 class trt2;
 model chg = y_t1 trt2  / solution covb;
 lsmeans trt2 / pdiff=all adjust=tukey;
 by _imputation_;
 ods output SolutionF=mixparms CovB=mixcovb LSMEANS=lsm01 DIFFS=lsmdiffs01;
run;
ods select all;

*# Summarize parameters;
proc mianalyze parms=mixparms;
 class trt2;
 modeleffects Intercept y_t1 trt2;
run;

proc means data=out01;
 class _imputation_ trt2;
 var y_t2;
run;

data _out01;
 set out01;
run;


*///////////////////////// STAGE 8 - MNAR test /////////////////////////;
*MNAR original;
proc mi data=lb05 out=out01 nimpute=500 seed=2012;
 class trt2 x2;
 var trt2 y_t1 x1 x2 y_t2;
 monotone reg(y_t2 / details); 
 mnar adjust(y_t2 / delta=-4 adjustobs=(trt2="A"));*deltas -1,-2,-3,-4,-5. After -4 starting to become significant;
run;

proc means data=out01;
 class _imputation_ trt2;
 var y_t2;
run;

data out02;
 set out01;
 chg = y_t2-y_t1;
run;

 proc sgplot data=out02;
  styleattrs datacontrastcolors=(blue red) datasymbols=(circleFilled);
  scatter x=y_t1 y=chg / group=trt2 ;
  where _imputation_=1;
 run;

ods select none;
proc mixed data=out02;
 class trt2;
 model chg = y_t1 trt2  / solution covb;
 lsmeans trt2 / pdiff=all adjust=tukey;
 by _imputation_;
 ods output SolutionF=mixparms CovB=mixcovb LSMEANS=lsm01 DIFFS=lsmdiffs01;
run;
ods select all;


*# Summarize parameters;
proc mianalyze parms=mixparms;
 class trt2;
 modeleffects Intercept y_t1 trt2;
run;

*# Summarizing Least Square Means;
proc mianalyze parms=lsm01;
 class trt2;
 modeleffects trt2;
 ods output parameterestimates=lsm02;
run;

*# Summarizing Least Square Means Differences;
proc mianalyze parms=lsmdiffs01;
 class trt2;
 modeleffects trt2;
 ods output parameterestimates=lsmdiffs02;
run;


*///////////////////////// STAGE 9 - modified MNAR test /////////////////////////;
*MNAR modified
* we try to create parameters for missings;
proc iml;

   nimpute= 500;
   call randseed( 2012);
   mean= { -5 -0.5};
   cov= { 0.5 0.01 , 0.01 0.05};

/*  ---- Simulate nimpute bivariate normal variates ---- */
   d= randnormal( nimpute, mean, cov);

   impu= j(nimpute, 1, 0);
   do j=1 to nimpute;  impu[j,]= j;  end;
   delta= impu || d;

/*  --- Output shift parameters for groups ---- */
   create parm1 from delta[colname={_Imputation_ Shift_P Shift_A}];
   append from delta;
quit;

proc mi data=lb05 out=out01 nimpute=500 seed=2012;
 class trt2 x2;
 var trt2 y_t1 x1 x2 y_t2;
 monotone reg(y_t2 / details); 
 mnar adjust(y_t2 / adjustobs=(trt2="A") parms(shift=shift_a)=parm1)
	  adjust(y_t2 / adjustobs=(trt2="P") parms(shift=shift_p)=parm1);
run;

proc means data=out01;
 class _imputation_ trt2;
 var y_t2;
run;

data out02;
 set out01;
 chg = y_t2-y_t1;
run;

 proc sgplot data=out02;
  styleattrs datacontrastcolors=(blue red) datasymbols=(circleFilled);
  scatter x=y_t1 y=chg / group=trt2 ;
  where _imputation_=1;
 run;

ods select none;
proc mixed data=out02;
 class trt2;
 model chg = y_t1 trt2  / solution covb;
 lsmeans trt2 / pdiff=all adjust=tukey;
 by _imputation_;
 ods output SolutionF=mixparms CovB=mixcovb LSMEANS=lsm01 DIFFS=lsmdiffs01;
run;
ods select all;


*# Summarize parameters;
proc mianalyze parms=mixparms;
 class trt2;
 modeleffects Intercept y_t1 trt2;
run;


*///////////////////////// STAGE 10 - analysis of trt factor /////////////////////////;
*means of change by trt group;
proc means data=fin_imput;
 by trt2;
 var chg;
run;

proc glm data=fin_imput;
class trt;
model chg = y_t1 trt;
run;