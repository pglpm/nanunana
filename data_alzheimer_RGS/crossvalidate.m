(* one-out cross-validation with eigendecomposition *)
ClearAll[crossvalidatet1e]
crossvalidatet1[data1_, otherdata_, pp_, setofquantities_] :=
  Block[{origsubjectfile,
	 numcategories = (Length@otherdata)+1,
	 numquantities,
	 allindividuals = Length[data1[[1]]],
	 numindividuals, tparameter, covmatrix, tmatrix, means,
	 results, quant, mainfile, secondfile,
	 distr, hits, totpr, datax},

	covmatrix[meansi_, datai_] := (datai - meansi).T[datai - meansi];

numindividuals[1]=allindividuals-1;

Do[origsubjectfile[ii] = otherdata[[ii-1]];
   numindividuals[ii] = Length[origsubjectfile[ii][[1]]];
   , {ii, 2, numcategories}];
	
ParallelTable[(* cycle through possible graph properties *)
  Block[{hits,tparameter,tmatrix,means,totpr=0},
	quant=setofquantities[[ll]];
	numquantities=Length[quant];

hits[1] = 0;
tparameter[1]=numindividuals[1]+1-numquantities;

Do[means[ii]=Mean@T[origsubjectfile[ii][[quant]]];
   tparameter[ii]=numindividuals[ii]+1-numquantities;
   tmatrix[ii]=covmatrix[means[ii], origsubjectfile[ii][[quant]]] *
			  (numindividuals[ii]+1)/numindividuals[ii]/tparameter[ii];
   hits[ii] = 0;
,{ii,2,numcategories}];
mainfile=data1[[quant]];

Do[(* cycle through individuals in first group *)
datax=T[mainfile][[ndat]];
secondfile=T[Drop[T[mainfile],{ndat}]];
means[1]=Mean@T[secondfile];
tmatrix[1]=covmatrix[means[1],secondfile] *
		      (numindividuals[1]+1)/numindividuals[1]/tparameter[1];

distr=Normalize[(* Normalize is faster *)
  Table[pp[[ii]]*
	PDF[MultivariateTDistribution[means[ii],tmatrix[ii],tparameter[ii]]
	   ,datax], {ii,numcategories}],Total];

++(hits[Last@Ordering[distr,-1]]);
totpr=totpr+distr;

,{ndat,allindividuals}];

{Array[hits,3],totpr}/allindividuals]

	, {ll, Length@setofquantities}, Method->"CoarsestGrained"]
];


(* one-out cross-validation *)
ClearAll[crossvalidatet1]
crossvalidatet1[data1_, otherdata_, pp_, setofquantities_] :=
  Block[{origsubjectfile,
	 numcategories = (Length@otherdata)+1,
	 numquantities,
	 allindividuals = Length[data1[[1]]],
	 numindividuals, tparameter, covmatrix, tmatrix, means,
	 results, quant, mainfile, secondfile,
	 distr, hits, totpr, datax},

	covmatrix[meansi_, datai_] := (datai - meansi).T[datai - meansi];

numindividuals[1]=allindividuals-1;

Do[origsubjectfile[ii] = otherdata[[ii-1]];
   numindividuals[ii] = Length[origsubjectfile[ii][[1]]];
   , {ii, 2, numcategories}];
	
ParallelTable[(* cycle through possible graph properties *)
  Block[{hits,tparameter,tmatrix,means,totpr=0},
	quant=setofquantities[[ll]];
	numquantities=Length[quant];

hits[1] = 0;
tparameter[1]=numindividuals[1]+1-numquantities;

Do[means[ii]=Mean@T[origsubjectfile[ii][[quant]]];
   tparameter[ii]=numindividuals[ii]+1-numquantities;
   tmatrix[ii]=covmatrix[means[ii], origsubjectfile[ii][[quant]]] *
			  (numindividuals[ii]+1)/numindividuals[ii]/tparameter[ii];
   hits[ii] = 0;
,{ii,2,numcategories}];
mainfile=data1[[quant]];

Do[(* cycle through individuals in first group *)
datax=T[mainfile][[ndat]];
secondfile=T[Drop[T[mainfile],{ndat}]];
means[1]=Mean@T[secondfile];
tmatrix[1]=covmatrix[means[1],secondfile] *
		      (numindividuals[1]+1)/numindividuals[1]/tparameter[1];

distr=Normalize[(* Normalize is faster *)
  Table[pp[[ii]]*
	PDF[MultivariateTDistribution[means[ii],tmatrix[ii],tparameter[ii]]
	   ,datax], {ii,numcategories}],Total];

++(hits[Last@Ordering[distr,-1]]);
totpr=totpr+distr;

,{ndat,allindividuals}];

{Array[hits,3],totpr}/allindividuals]

	, {ll, Length@setofquantities}, Method->"CoarsestGrained"]
];


ClearAll[crossvalidatetall];
crossvalidatetall[dataall_, pp_, setofquantities_] :=
  Block[{numcategories = Length@dataall},
Table[crossvalidatet1[dataall[[ii]], Drop[dataall,{ii}], pp, setofquantities], {ii,ncat}]];





(* not efficient *)
ClearAll[crossvalidate1old];
crossvalidate1old[data1_, otherdata_, pp_, prior_] :=
  Block[{ncat = Length@otherdata+1,
	 nqua = Length@data1,
	 nind1 = Length[data1[[1]]] - 1,
	 nind, means, tmatrix, scmatrix, data2, tparameter, datum,
	 onind, omeans, distr, hits = 0, totpr = 0,
	 rnind, rmeans, rscmatr
    },

    scmatrix[meansi_, datai_] := (datai - meansi).T[datai - meansi];

If[prior == 0,
   rnind = 0; rmeans = Table[0, {nqua}]; rscmatr = Table[0, {nqua}, {nqua}];,
 {rnind, rmeans, rscmatr} = prior; rscmatr = rscmatr * rnind;
];

Do[(* combine prior data + reference prior *)
  onind = Length[otherdata[[ii-1]][[1]]];
  omeans = Mean@T[otherdata[[ii-1]]];
  nind[ii] = onind + rnind;
  tparameter[ii] = nind[ii] + 1 - nqua;
  means[ii] = (rnind * rmeans + onind * omeans)/nind[ii];
  tmatrix[ii] = Symmetrize[scmatrix[omeans, otherdata[[ii-1]]] + rscmatr +
			   T[{omeans - rmeans}].{omeans-rmeans} *
			   onind * rnind/nind[ii]] *
		 (nind[ii]+1)/(nind[ii]*tparameter[ii]);
, {ii, 2, ncat}
];

nind[1] = nind1 + rnind;
tparameter[1] = nind[1] + 1 - nqua;

    If[prior==0,
Do[(* cycle through individuals in first group *)
  datum = T[data1][[idat]];

  data2 = T[Drop[T[data1], {idat}]];
  omeans = Mean@T[data2];
   means[1] = omeans;
  tmatrix[1] = Symmetrize[scmatrix[omeans, data2] *
		 (nind[1]+1)/(nind[1]*tparameter[1])];
  
  distr = Normalize[(* Normalize is faster *)
    Table[pp[[ii]] *
          PDF[
            MultivariateTDistribution[means[ii], tmatrix[ii], tparameter[ii]],
	    datum]
	, {ii, ncat}], Total];
     
     If[Last@Ordering[distr, -1] == 1, ++hits, Null];
     totpr = totpr + distr[[1]];
     
     , {idat, nind1 + 1}];
       ,
    Do[(* cycle through individuals in first group *)
  datum = T[data1][[idat]];

  data2 = T[Drop[T[data1], {idat}]];
  omeans = Mean@T[data2];
  means[1] = (rnind * rmeans + nind1 * omeans)/nind[1];
  tmatrix[1] = Symmetrize[scmatrix[omeans, data2] + rscmatr +
			   T[{omeans - rmeans}].{omeans-rmeans} *
			   nind1 * rnind/nind[1]] *
		 (nind[1]+1)/(nind[1]*tparameter[1]);
  
  distr = Normalize[(* Normalize is faster *)
    Table[pp[[ii]] *
          PDF[
            MultivariateTDistribution[means[ii], tmatrix[ii], tparameter[ii]],
	    datum]
	, {ii, ncat}], Total];
     
     If[Last@Ordering[distr, -1] == 1, ++hits, Null];
     totpr = totpr + distr[[1]];
     
     , {idat, nind1 + 1}];
]

    
    {hits, totpr}/(nind1 + 1)
 ];

ClearAll[crossvalidateallold];
crossvalidateall[dataall_, pp_, prior_] :=
  Block[{ncat = Length@dataall},
Table[crossvalidate1old[dataall[[ii]], Drop[dataall,{ii}], pp, prior], {ii,ncat}]
  ];
