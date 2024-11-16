*Step 1: Importing the files to library;
PROC IMPORT DATAFILE='Duration\Project\liver_cirrhosis.csv'
    OUT=liver_cirrhosis
    DBMS=CSV
    REPLACE;
    GETNAMES=YES;
RUN;

PROC CONTENTS DATA=liver_cirrhosis;
RUN;

PROC PRINT DATA=liver_cirrhosis (OBS=10);
RUN;

*data has no missing values. The data is complete and the data
imputation is not needed;
*here we remove rows where the censoring due to transplantation is made;
Proc SQL;
    Create table LC_2 as
    select *,
    case when Status = "D" then 1 else 0 end AS StatusNum,
    case when Age > 14600 THEN 'Older than 40' ELSE '40 or Younger' end
as Age_group
    from liver_cirrhosis
    where Status in ("C","D");
Quit;

*Descriptive analysis;
proc means data=LC_2;
/*class statusnum;*/
run;

*Because in the data description it is said that Status has three
distinct values, we check for it.
As it turns out, there are in fact only two values: C and D;
proc freq data=LC_2;
TABLES Status Sex Drug;
run;

*some more descriptive histograms to understand distribution of data;
PROC SGPLOT DATA=LC_2;
    HISTOGRAM Age;
    DENSITY Age;
    TITLE "Age Distribution";
RUN;

PROC SGPLOT DATA=LC_2;
    VBAR Sex / GROUP=Status;
    TITLE "Survival Status by Sex";
RUN;

*Step 3: Non-Parametric Survival Analysis (Kaplan-Meier Estimator);
PROC LIFETEST DATA=LC_2 method=lt PLOTS=(SURVIVAL, s, h, p); *notable,
for alternative -- put this and remove method=lt;
    TIME N_Days*StatusNum(0);
RUN;

*categorical variable Drug;
PROC LIFETEST DATA=LC_2 method=lt PLOTS=(SURVIVAL, s, h, p); *notable,
for alternative -- put this and remove method=lt;
    TIME N_Days*StatusNum(0);
    STRATA Drug;
RUN;

*test of hypothesis 1;
PROC LIFETEST DATA=LC_2 method=lt PLOTS=(SURVIVAL, s, h, p); *notable,
for alternative -- put this and remove method=lt;
    TIME N_Days*StatusNum(0);
    STRATA age_group;
RUN;

PROC LIFETEST DATA=LC_2 method=lt PLOTS=(SURVIVAL, s, h, p); *notable,
for alternative -- put this and remove method=lt;
    TIME N_Days*StatusNum(0);
    STRATA age_group Drug;
RUN;

PROC LIFETEST DATA=LC_2 method=lt PLOTS=(SURVIVAL, s, h, p); *notable,
for alternative -- put this and remove method=lt;
    TIME N_Days*StatusNum(0);
    STRATA Sex;
RUN;

*Step 3.2 NP Kaplan-Meier;
PROC LIFETEST DATA=LC_2 method=pl notable PLOTS=(s, ls, lls);
    TIME N_Days*StatusNum(0);
    /*STRATA Drug Age_group Sex;*/
RUN;

*NP Kaplan-Meier with categorical variables;
PROC LIFETEST DATA=LC_2 method=pl notable PLOTS=(s, ls, lls);
    TIME N_Days*StatusNum(0);
    STRATA Drug;
RUN;

*NP Kaplan-Meier with categorical variables;
PROC LIFETEST DATA=LC_2 method=pl notable PLOTS=(s, ls, lls);
    TIME N_Days*StatusNum(0);
    STRATA Age_group;
RUN;

*NP Kaplan-Meier with categorical variables;
PROC LIFETEST DATA=LC_2 method=pl notable PLOTS=(s, ls, lls);
    TIME N_Days*StatusNum(0);
    STRATA Drug Age_group;
RUN;

*NP Kaplan-Meier with categorical variables;
PROC LIFETEST DATA=LC_2 method=pl notable PLOTS=(s, ls, lls);
    TIME N_Days*StatusNum(0);
    STRATA Sex;
RUN;

*to add all variables needed for testing our hypotheses with NP Models;
PROC LIFETEST DATA=LC_2 method=pl notable PLOTS=(s, ls, lls); *notable,
for alternative -- put this and remove method=lt;
    TIME N_Days*StatusNum(0);
    STRATA Drug Age_group Sex;
RUN;

*Part4 Parametric Models.
Classical approach. Exponential model
Part 4.1.;
PROC LIFEREG DATA=LC_2;
    CLASS Drug Age_group Sex;
    MODEL N_Days*StatusNum(0) = / DISTRIBUTION=Exponential;
RUN;

PROC LIFEREG DATA=LC_2;
    CLASS Drug Age_group Sex;
    MODEL N_Days*StatusNum(0) = Age Age_group Sex Bilirubin Albumin
Prothrombin Drug / DISTRIBUTION=Exponential;
RUN;

*Part4.2. Weibull Model
without explanatory;
PROC LIFEREG DATA=LC_2;
    CLASS Drug Age_group Sex;
    MODEL N_Days*StatusNum(0) = / DISTRIBUTION=Weibull;
RUN;

*with explanatory;
PROC LIFEREG DATA=LC_2;
    CLASS Drug Age_group Sex;
    MODEL N_Days*StatusNum(0) = Age Age_group Sex Bilirubin Albumin
Prothrombin Drug / DISTRIBUTION=Weibull;
RUN;

*Part4.3. Log-Normal Model;
PROC LIFEREG DATA=LC_2;
    CLASS Drug Age_group Sex;
    MODEL N_Days*StatusNum(0) = / DISTRIBUTION=Lognormal;
RUN;

PROC LIFEREG DATA=LC_2;
    CLASS Drug Age_group Sex;
    MODEL N_Days*StatusNum(0) = Age Age_group Sex Bilirubin Albumin
Prothrombin Drug / DISTRIBUTION=Lognormal;
RUN;

*Part5 Parametric Models Bayesias approach
Part5.1. Without explanatory variables
Weibul model;
PROC LIFEREG DATA=LC_2;
    CLASS Drug Age_group Sex;
    MODEL N_Days*StatusNum(0) = / DISTRIBUTION=Weibull;
    bayes;
    ods output PosteriorSample=PS_WEIBULL;
RUN;

*Log-Normal Model;
PROC LIFEREG DATA=LC_2;
    CLASS Drug Age_group Sex;
    MODEL N_Days*StatusNum(0) = / DISTRIBUTION=LNORMAL;
    bayes;
    ods output PosteriorSample=PS_LNORMAL;
RUN;

*Part 5.2. With explanatory variables;
PROC LIFEREG DATA=LC_2;
    CLASS Sex;
    MODEL N_Days*StatusNum(0) = Age Sex / DISTRIBUTION=Weibull;
    bayes nbi=2000 nmc=10000 THIN=2 coeffprior=normal diagnostics=all;
    ods output PosteriorSample=PS_WEIBULL;
RUN;

*Log-Normal Model;
PROC LIFEREG DATA=LC_2;
    CLASS Sex;
    MODEL N_Days*StatusNum(0) = Age Sex / DISTRIBUTION=LNORMAL;
    bayes nbi=2000 nmc=10000 THIN=2 coeffprior=normal diagnostics=all;
    ods output PosteriorSample=PS_LNORMAL;
RUN;