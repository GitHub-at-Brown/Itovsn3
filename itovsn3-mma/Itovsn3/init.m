(* ::Package:: *)

(*
    Itovsn3.m: a Mathematica package for
    Symbolic Ito calculus Itovsn3. Version 3.70,
      Copyright March, June 1992, April 1998, November 2002, January 2003
    Author:   Wilfrid S.Kendall.

    Version 3.68 was tested on 
    Mathematica version 1, Macintosh, Windows 3.1;
    Mathematica version 2, Windows 3.1;
    Mathematica version 2.2.3, Windows 97;
    Mathematica version 3.0, Linux.
  
    Version 3.69 has been tested on
    Mathematica version 4.0.0.0, Windows 2000

    Version 3.70 has been tested on
    Mathematica version 4.0.0.0, Windows 2000
    Mathematica version 4.1.0.0, Unix

The author disclaims any responsibility for problems experienced by
users as a result of using this software, but welcomes your feedback!
The software may not be copied or distributed for commercial gain
without written permission from the author, but may be passed to
others if left intact in its present form.

Comments please to  w.s.kendall@warwick.ac.uk  or to

	Prof W S Kendall,
	Department of Statistics,
	University of Warwick,
	Coventry CV4 7AL, UK.
*)

BeginPackage["Itovsn3`"]

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
"RandomQ[x,sdl]==True  if  x  is an expression in semimartingales or\n
stochastic differentials excluding those given in  sdl."

BSDQ::usage =
"BSDQ[x,sdl]==True  if  x  is a basic stochastic differential\n
excluding those given in  sdl."

GetItoProc::usage =
"GetItoProc[]  gives a list of all defined semimartingales."
(* ============================= 					*)

Begin["`private`"]

(* =>	Implementing the Ito formula:					*)
ItoD[f_] := Block[
		{ ff =(f/.ItoIntegral[sdx_]->ItoIntegral[sdx,t]) },
		ItoExpand[ 
                 (Dt[ff,t] dt + (1/2) Dt[ff,{t,2}] dt^2)/.Freeze[sd_]->sd ] 
		];
ItoExpand[sd_] := ((Expand[sd] 
			/. ItoIntegral[sdx_,t]->ItoIntegral[sdx])
			/. ItoMultiplications);

(* =>	Implementing Drift:						*)
Drift[sd_] := Apply[Plus,Map[Coefficient[Expand[sd],#] Drftbydt[#] dt &,CSD]];

(* =>	Updating first- and second-order structure:			*)
AddQuadVar[dX_ dY_,xdt_] := AddQuadVar[dX,dY,xdt];
AddQuadVar[dX_^2,xdt_]   := AddQuadVar[dX,dX,xdt];
AddQuadVar[dX_,dY_,xdt_] := (Brktbydt[dX,dY] = xdt/dt);
AddDrift[dX_,DriftdX_]   := (Drftbydt[dX] = DriftdX/dt);

(* =>	Finding initial value:						*)
AddFixed[t0_,y_,y0_] := (Fixed[t0,y]=y0);
Derivative[1,0][Fixed][t0_,x_] := 0;
Derivative[0,1][Fixed][t0_,x_] := 0;
InitialValue[t0_,x_] :=
	((((x /. Map[#->MayFix[t0,#]&,CS])
	/. ItoIntegral[y_] ->MayFix[t0,ItoIntegral[y]])
	/. MayFix[a_,y_] -> MayFix[a,(y/.MayFix[t0,u_]->u)])
	/. MayFix[t0,u_] -> Fixed[t0,u]);

(* =>	Simplest properties of ItoIntegral:				*)
Derivative[0,2][ItoIntegral][sd_,t_] := 0;
Derivative[0,1][ItoIntegral][sd_,t_] := Freeze[sd]/ItoD[t];
Derivative[1][Freeze][sd_] := 0;
Derivative[1,0][ItoIntegral][sd_,t_] := 0;

(* =>	Introducing basic semimartingale
	and associated stochastic differential:				*)
Introduce[X_, dX_] :=
	If[Not[MemberQ[CS,X]] && Not[MemberQ[CSD,dX]],
	(
	X/:  Dt[X,t] = dX/dt;
	dX/: Dt[dX,t] = 0 ;
        ItoIntegral[dX] = X - Fixed[0,X];
	CSD = Prepend[CSD,dX];
	CS = Prepend[CS,X];
	ItoMultiplications = Join[
		Map[(dX # -> Brktbydt[dX,#] dt)&,CSD],
			ItoMultiplications	];
	AddQuadVar[dX dt,0];
	CSD
	), "Attempt to re-introduce semimartingale or stochastic differential!"
	];

(* =>	Initialization:							*)
ItoInit[time_,dtime_] :=
	(
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

(* =>	Resetting stochastic calculus structures:			*)
ItoReset[time_,dtime_] :=
	(
	Map[ItoClear[#,ItoD[#]]&,CS];
	Clear[Fixed];
	Clear[Brktbydt];
	Clear[Drftbydt];
        Clear[ItoIntegral];
	TableForm[Prepend[First[ItoInit[time,dtime]],"Itovsn3  resetting ..."]]
	);

(* =>	Test if x involves semimartingales or differentials except 
        those given in the second list (presumed deterministic!)*)
RandomQ[x_,sdl_] := 
        Not[FreeQ[x,Apply[Alternatives,Complement[Union[CS,CSD],sdl]]]]

(* =>	Test if x is a basic stochastic differential excepting
        those given in the second list*)
BSDQ[x_,sdl_] := 
        MemberQ[Complement[CSD,sdl],x]

(* =>	Clear aspects of semimartingale  X  and differential  dX:	*)
ItoClear[X_,dX_] :=
	( dX/: Dt[dX,t] = .; X/:  Dt[X,t] = .; );

(* =>	Reporting status:						*)
ItoStatus[] :=
	Print[ColumnForm[{
        "---------------",
	"Summary of current structure of stochastic differentials",
	"- - - - - - - -",
	"Current second-order structure of semimartingale differentials:",
	TableForm[Outer[ItoExpand[#1 #2]&,CSD,CSD],TableHeadings->{CSD,CSD}],
	"- - - - - - - -",
	"Current first-order structure of semimartingale differentials:",
	TableForm[Map[Drift[#]&,CSD],
		  TableDirections->Row,TableHeadings->{CSD,{"Drifts:"}}],
	"- - - - - - - -",
	"Current initial values:",
	TableForm[Map[Fixed[0,#]&,CS],
		TableDirections->Row,TableHeadings->{CS,{"Initially:"}}],
	"---------------"
	}]];

(* =>	Brownian Basis:							*)
BrownBasis[SemimartingaleList_,InitialValueList_] :=
	(
	Map[Apply[BrownSingle,#]&,
	Transpose[{SemimartingaleList,InitialValueList}]];
	BrownPairs[Map[ItoD,SemimartingaleList]]
	);
BrownSingle[X_,X0_] :=
	Block[
		{dX=ToExpression[StringJoin["d",ToString[X]]]},
		Introduce[X,dX];
		AddQuadVar[dX^2,dt];
		AddDrift[dX,0];
		Fixed[0,X]=X0
	];
BrownPairs[SL_] :=
	If[Length[SL]>1,
	Map[AddQuadVar[First[SL] #,0]&,Rest[SL]];BrownPairs[Rest[SL]]];

(* =>   Definition using Ito stochastic differential equations:		*)
Itosde[X_,dX_==sd_,X0_] :=
	(
	Introduce[X,dX];
	AddDrift[dX,Drift[sd]];
	Fixed[0,X]=X0;
	AddQuadVar[dX^2,ItoExpand[sd^2]];
	Map[AddQuadVar[dX #,ItoExpand[sd #]]&,CSD];
	);     		

(* => List defined semimartingales and differentials		*)
GetItoProc[]:={CS,CSD};

End[]
(* =============================					 *)

EndPackage[]
(* =============================					 *)
