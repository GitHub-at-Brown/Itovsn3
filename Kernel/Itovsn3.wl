(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["FernandoDuarte`Itovsn3`"];


(* ::Subsection:: *)
(*Public symbols*)


AddDrift::usage =
"AddDrift[dX,DriftdX]  sets  Drftbydt[dX] = DriftdX/dt,  so that\n
(for example)  Drift[dX]  yields  DriftdX."

AddFixed::usage =
"AddFixed[t0,X,X0]  sets  Fixed[t0,X] = X0,  so that\n
(for example)  InitialValue[t0,X]  yields  X0."

AddQuadVar::usage =
"AddQuadVar[dX,dY,xdt]  sets  Brktbydt[dX,dY] = xdt/dt,  altering \n
the effect of substitution using  ItoMultiplications,  so that\n
(for example)  ItoExpand[dX dY]  yields  xdt.  Similarly for\n
AddQuadVar[dX dY,xdt]  and  AddQuadVar[dX^2,xdt]."

Brktbydt::usage =
"Brktbydt[dX,dY]  is a placeholder for the FORMAL quotient by  dt  of\n
the bracket differential  dX dY."

BrownBasis::usage =
"BrownBasis[SemimartingaleList,InitialValueList]  introduces and sets\n
up second- and first-order structure for an independent set of Brownian\n
basic semimartingales, identifiers in  SemimartingaleList,  with initial \n
value expressions given by  InitialValueList.  It uses  BrownSingle[X,X0] \n
as a supplementary procedure. Corresponding basic stochastic differential\n
identifiers are created by prepending 'd' to the  SemimartingaleList  \n
identifiers."

BrownSingle::usage =
"BrownSingle[X,X0]  introduces and sets up a single Brownian basic\n
semimartingale identifier  X  (creating the semimartingale differential\n
identifier  dX)  with initial value expression  X0."
 
Drftbydt::usage =
"Drftbydt[dX]  is a placeholder for the FORMAL quotient by  dt
of Drift of  dX.\n"

Drift::usage =
"Drift[sd]  computes the drift differential of the stochastic\n
differential expression  sd.  NOTE that  Drift  assumes that  sd  is\n
genuinely a stochastic differential expression."

Fixed::usage =
"Fixed[t0,X]  is a placeholder for the fixed value of the basic \n
semimartingale X  at time  t0  (usually  t0=0)."

InitialValue::usage =
"InitialValue[t0,f]  computes the value of the expression  f  \n
at time  t0."

Introduce::usage =
"Introduce[Smgl,dSmgl]  introduces the basic semimartingale identifier \n
Smgl  with associated basic stochastic differential identifer  dSmgl.  \n
Attempts to reintroduce semimartingale or stochastic differential\n
identifiers are reported as errors."

ItoD::usage =
"ItoD[f]  computes the stochastic differential of the semimartingale\n
expression  f.  NOTE that  ItoD  assumes that  f  is a genuine\n
semimartingale expression!"

ItoExpand::usage =
"ItoExpand[sd]  computes simplication of the stochastic differential\n
products in the expression  sd  which is its argument."

ItoInit::usage =
"ItoInit[t,dt]  starts things off with basic structures, using the \n
identifier  t  for time variable and the identifier  dt  for its \n
differential."

ItoIntegral::usage =
"ItoIntegral[sd]  represents the Ito integral of the stochastic\n
differential expression  sd."

ItoReset::usage =
"ItoReset[t,dt]  resets all structures, using  ItoInit[t,dt]."

Itosde::usage =
"Itosde[X,dX==sd,X0]  introduces and sets up a basic semimartingale \n
identifier X  with basic stochastic differential identifier  dX  and \n
initial value expression  X0,  and satisfying the second- and first- \n
order structure implied by the stochastic differential equation  dX==sd."

ItoStatus::usage =
"ItoStatus[]  reports current structures."

RandomQ::usage =
"RandomQ[x]==True if x is an expression in semimartingales or stochastic differentials."<>"\n"<>
"RandomQ[x,sdl]==True  if  x  is an expression in semimartingales or\n
stochastic differentials excluding those given in  sdl."

BSDQ::usage =
"BSDQ[x,sdl]==True  if  x  is a basic stochastic differential\n
excluding those given in  sdl."

GetItoProc::usage =
"GetItoProc[]  gives a list of all defined semimartingales."


(* ::Section:: *)
(*Code*)


Begin["`Private`"];


(* ::Subsubsection:: *)
(*Implementing  the  Ito  formula*)


ItoD[f_] := Block[
		{ ff =(f/.ItoIntegral[sdx_]:>ItoIntegral[sdx,t]) },
		ItoExpand[ 
                 (Dt[ff,t] dt + (1/2) Dt[ff,{t,2}] dt^2)/.Freeze[sd_]:>sd ] 
		];
ItoExpand[sd_] := ((Expand[sd] 
			/. ItoIntegral[sdx_,t]:>ItoIntegral[sdx])
			/. ItoMultiplications);


(* ::Subsubsection:: *)
(*Implementing  Drift*)


Drift[sd_] := Apply[Plus,Map[Coefficient[Expand[sd],#] Drftbydt[#] dt &,CSD]];


(* ::Subsubsection:: *)
(*Updating  first - and  second - order  structure*)


AddQuadVar[dX_ * dY_,xdt_] := AddQuadVar[dX,dY,xdt];
AddQuadVar[dX_^2,xdt_]   := AddQuadVar[dX,dX,xdt];
AddQuadVar[dX_,dY_,xdt_] := (Brktbydt[dX,dY] = xdt/dt);
AddDrift[dX_,driftdX_]   := (Drftbydt[dX] = driftdX/dt);


(* ::Subsubsection:: *)
(*Finding  initial  value*)


AddFixed[t0_,y_,y0_] := (Fixed[t0,y]=y0);
Derivative[1,0][Fixed][t0_,x_] := 0;
Derivative[0,1][Fixed][t0_,x_] := 0;
InitialValue[t0_,x_] :=
	((((x /. Map[#->MayFix[t0,#]&,CS])
	/. ItoIntegral[y_] :> MayFix[t0,ItoIntegral[y]])
	/. MayFix[a_,y_] :> MayFix[a,(y/.MayFix[t0,u_]:>u)])
	/. MayFix[t0,u_] :> Fixed[t0,u]);


(* ::Subsubsection:: *)
(*Simplest  properties  of  ItoIntegral*)


Derivative[0,2][ItoIntegral][sd_,t_] := 0;
Derivative[0,1][ItoIntegral][sd_,t_] := Freeze[sd]/ItoD[t];
Derivative[1][Freeze][sd_] := 0;
Derivative[1,0][ItoIntegral][sd_,t_] := 0;


(* ::Subsubsection:: *)
(*Introducing  basic  semimartingale  and  associated  stochastic  differential*)


Introduce[x_, dx_] := If[
	Not[MemberQ[CS,x]] && Not[MemberQ[CSD,dx]],
	(
		x/:  Dt[x,t] = dx/dt;
		dx/: Dt[dx,t] = 0 ;
	        ItoIntegral[dx] = x - Fixed[0,x];
		CSD = Prepend[CSD,dx];
		CS = Prepend[CS,x];
		ItoMultiplications = Join[
			Map[(dx # -> Brktbydt[dx,#] dt)&,CSD],
				ItoMultiplications	];
		AddQuadVar[dx dt,0];
		CSD
	),
	"Attempt to re-introduce semimartingale or stochastic differential."
];


(* ::Subsubsection:: *)
(*Initialization*)


ItoInit[time_,dtime_] :=(
	t = time;
	dt = dtime;
	SetAttributes[Brktbydt,Orderless];
	CS = {};
	CSD = {};
	ItoMultiplications = {};
	Introduce[t,dt];
	AddDrift[dt,dt];
	ItoIntegral[0] = 0; 
	Fixed[t0_,t] = t0;
	Fixed[0,ItoIntegral[y_]] = 0;
	TableForm[{
		"Itovsn3  initialized",
		SequenceForm["with time semimartingale ",t],
		SequenceForm["and time differential ",dt]
		}] 
);


(* ::Subsubsection:: *)
(*Resetting  stochastic  calculus  structures*)


ItoReset[time_,dtime_] :=(
	Map[ItoClear[#,ItoD[#]]&,CS];
	Clear[Fixed];
	Clear[Brktbydt];
	Clear[Drftbydt];
        Clear[ItoIntegral];
	TableForm[Prepend[First[ItoInit[time,dtime]],"Itovsn3  resetting ..."]]
);


(* ::Subsubsection:: *)
(*Test  if  x  involves  semimartingales  or  differentials, except those  given  in  the  second  list (if any)*)


RandomQ[x_,sdl_:{}] := Not[FreeQ[x,Apply[Alternatives,Complement[Union[CS,CSD],sdl]]]]


(* ::Subsubsection:: *)
(*Test  if  x  is  a  basic  stochastic  differential  excepting  those  given  in  the  second  list  (if any)*)


BSDQ[x_,sdl_:{}] := MemberQ[Complement[CSD,sdl],x]


(* ::Subsubsection:: *)
(*Clear  aspects  of  semimartingale   x   and  differential   dx*)


ItoClear[x_,dx_] := ( dx/: Dt[dx,t] =. ; x/:  Dt[x,t] =. ; );


(* ::Subsubsection:: *)
(*Reporting  status*)


ItoStatus[opts : OptionsPattern[{Grid}]] := Grid[{
	{" Summary of current structure of stochastic differentials"},
	{" Current second-order structure of semimartingale differentials:"},
	{
		TableForm[Outer[ItoExpand[#1 #2]&,CSD,CSD],TableHeadings->{CSD,CSD}]
	},
	{" Current first-order structure of semimartingale differentials:"},
	{
		TableForm[Map[Drift[#]&,CSD],
		  TableDirections->Row,TableHeadings->{CSD,{"Drifts:"}}]
	},
	{" Current initial values:"},
	{
		TableForm[Map[Fixed[0,#]&,CS],
		TableDirections->Row,TableHeadings->{CS,{"Initially:"}}]
	}
	},
	opts,
	Dividers->{{False},{True,True,False,False,False,False,False,True}},
	Spacings->{Automatic,1.5},
	Background->{None,{1->LightGray,2->LightBlue,3->White,4->LightBlue,5->White,6->LightBlue}}
]


(* ::Subsubsection:: *)
(*Brownian  Basis*)


BrownBasis[ semimartingaleList_,initialValueList_] :=(
	Map[Apply[BrownSingle,#]&,
	Transpose[{ semimartingaleList,initialValueList}]];
	BrownPairs[Map[ItoD, semimartingaleList]]
);
BrownSingle[x_,x0_] :=
	Block[
		{dx=ToExpression[StringJoin["d",ToString[x]]]},
		Introduce[x,dx];
		AddQuadVar[dx^2,dt];
		AddDrift[dx,0];
		Fixed[0,x]=x0
	];
BrownPairs[sL_] :=
	If[Length[sL]>1,
	Map[AddQuadVar[First[sL] #,0]&,Rest[sL]];BrownPairs[Rest[sL]]];


(* ::Subsubsection:: *)
(*Definition  using  Ito  stochastic  differential  equations*)


Itosde[x_,dx_==sd_,x0_] :=(
	Introduce[x,dx];
	AddDrift[dx,Drift[sd]];
	Fixed[0,x]=x0;
	AddQuadVar[dx^2,ItoExpand[sd^2]];
	Map[AddQuadVar[dx #,ItoExpand[sd #]]&,CSD];
);     		


(* ::Subsubsection:: *)
(*List  defined  semimartingales  and  differentials*)


GetItoProc[]:={CS,CSD};


(* ::Section:: *)
(*Package Footer*)


End[];
EndPackage[];
