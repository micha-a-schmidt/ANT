(* ::Package:: *)

Print["ANT Version: 1.2 ( 08 Nov 2015)."];
Print["Authors: Paul W. Angel, Yi Cai, Nicholas L. Rodd, Michael A. Schmidt, Raymond R. Volkas"];
Print["Please cite: arXiv:1308.0463."];
Print["http://ant.hepforge.org"];
Print["The package ANT is written for Mathematica 8 and higher and it is distributed under the terms of GNU Public License http://www.gnu.org/copyleft/gpl.html"];
Print["Use ANT[expr] to evaluate the Passarino-Veltman functions contained in expr to its finite part in the limit of vanishing external momenta. Besides the function ANT, it defines the functions A0ant, B0ant, C0ant as well as D0ant, which take the same arguments as A0, B0i, C0i and D0i, respectively, and return the finite part of the PV function in the limit of external momenta."];
Print["Note that the package ANT has to be loaded after FormCalc"];


BeginPackage["ANT`"];

{ANT,A0ant,B0ant,C0ant,D0ant, AFant,AFDant,BFant,BFDant,CFant,CFDant,DFant,DFDant}

ANT::usage="The function ANT can be applied to any expression expr, e.g. ANT[expr], and it replaces all Passarino-Veltman functions (as defined by FormCalc/LoopTools) by the evaluated expression in the limit of vanishing external momenta. For the PV functions B0, B1, C0, C1, C2, C00, it also returns the expansion of the PV function up to the first derivatives with respect to the momenta."
A0ant::usage="A0ant[\!\(\*SuperscriptBox[\(m\), \(2\)]\)] returns the finite part of \!\(\*SubscriptBox[\(A\), \(0\)]\)(\!\(\*SuperscriptBox[\(m\), \(2\)]\))";
B0ant::usage="B0ant[id,\!\(\*SuperscriptBox[\(p\), \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(1\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(2\)], \(2\)]\)] returns the finite part of \!\(\*SubscriptBox[\(B\), \(id\)]\)(\!\(\*SuperscriptBox[\(p\), \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(1\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(2\)], \(2\)]\)) in the limit of vanishing external momenta. Additionally, it returns the first derivative of \!\(\*SubscriptBox[\(B\), \(0\)]\) as well as \!\(\*SubscriptBox[\(B\), \(1\)]\).";
C0ant::usage="C0ant[id,\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(1\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(2\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(12\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(1\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(2\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(3\)], \(2\)]\)] returns the finite part of \!\(\*SubscriptBox[\(C\), \(id\)]\)(\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(1\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(2\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(12\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(1\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(2\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(3\)], \(2\)]\)) in the limit of vanishing external momenta. Additionally, it returns the first derivative of \!\(\*SubscriptBox[\(C\), \(0\)]\), \!\(\*SubscriptBox[\(C\), \(1\)]\), \!\(\*SubscriptBox[\(C\), \(2\)]\) as well as \!\(\*SubscriptBox[\(C\), \(00\)]\).";
D0ant::usage="D0ant[id,\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(1\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(2\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(3\)], \(2\)]\), \!\(\*SuperscriptBox[SubscriptBox[\(p\), \(4\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(12\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(23\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(1\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(2\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(3\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(4\)], \(2\)]\)] returns the finite part of \!\(\*SubscriptBox[\(D\), \(id\)]\)(\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(1\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(2\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(3\)], \(2\)]\), \!\(\*SuperscriptBox[SubscriptBox[\(p\), \(4\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(12\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(23\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(1\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(2\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(3\)], \(2\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(m\), \(4\)], \(2\)]\)) in the limit of vanishing external momenta.";

EndPackage[];



(* ::Section:: *)
(*Definitions*)


(* ::Subsection::Closed:: *)
(*A*)


AFant[ms_]:=ms(1-Log[ms]);
AFant[0]:=0;
AFDant[ms_]:=0;


(* ::Subsection::Closed:: *)
(*B*)


BFant[LoopTools`bb0,m1s_,m2s_]:=
If[SameQ[m1s,m2s],-Log[m1s],1-(m1s Log[m1s]-m2s Log[m2s])/(m1s-m2s)];
BFDant[LoopTools`bb0,ps_,m1s_,m2s_]:=Block[{t},
t=m1s/m2s;If[SameQ[m1s,m2s],ps/(6 m2s),ps/m2s (-1+t^2-2 t Log[t])/(2 (-1+t)^3)]];


BFant[LoopTools`bb1,m1s_,m2s_]:=Block[{t},
t=m2s/m1s;If[SameQ[m1s,m2s],Log[m1s]/2,1/2 Log[m1s]+(-3+4t-t^2-4t Log[t]+2t^2 Log[t])/(4(-1+t)^2)]];
BFDant[LoopTools`bb1,ps_,m1s_,m2s_]:=Block[{t},
t=m1s/m2s;If[SameQ[m1s,m2s],-(ps/(12m2s)),ps/m2s (-(-1+t) (-1+t (5+2 t))+6 t^2 Log[t])/(6 (-1+t)^4)]];


BFant[LoopTools`bb00,m1s_,m2s_]:=Block[{t},
t=m2s/m1s;
If[SameQ[m1s,m2s],-(1/2)m1s(-1+Log[m1s]),1/8 m1s (3+3 t-2 (1+t) Log[m1s]-(2 t^2 Log[t])/(-1+t))]];




BFant[LoopTools`bb11,m1s_,m2s_]:=Block[{t},t=m2s/m1s;
If[m1s===m2s,-(Log[m1s]/3),-(1/3)Log[m1s]-(-(-1+t) (11+t (-7+2 t))+6 t (3+(-3+t) t) Log[t])/(18 (-1+t)^3)]];


BFDant[id_,ps_,m1s_,m2s_]:=Block[{},
Print["Derivative for B0i[",id,"] not defined"];
B0i[id,ps,m1s,m2s]-BFant[id,m1s,m2s]
];


(* ::Subsection:: *)
(*C*)


CFant[LoopTools`cc0,m1s_,m2s_,m3s_]:=Block[{t1,t2,s123,s12,s13,s23,us},
t1=m1s/m3s;t2=m2s/m3s;
s123=SameQ[m1s,m2s,m3s];
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s];
s13=SameQ[m1s,m3s]&&UnsameQ[m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m1s,m3s];
us=UnsameQ[m1s,m2s,m3s];
Which[
s123, -(1/(2m3s)),
s12,(1-t1+Log[t1])/(m3s (-1+t1)^2),
s13,(-1+t2-t2 Log[t2])/(m3s (-1+t2)^2),
s23,(-1+t1-t1 Log[t1])/(m3s (-1+t1)^2),
us, 1/m3s ((t1-t1 t2) Log[t1]+(-1+t1) t2 Log[t2])/((-1+t1) (t1-t2) (-1+t2))]
];
CFDant[LoopTools`cc0,p1s_,p2s_,p12s_,m1s_,m2s_,m3s_]:=Block[{t1,t2,s123,s12,s13,s23,us},
t1=m1s/m3s;t2=m2s/m3s;
s123=SameQ[m1s,m2s,m3s];
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s];
s13=SameQ[m1s,m3s]&&UnsameQ[m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m1s,m3s];
us=UnsameQ[m1s,m2s,m3s];
Which[
s123,-((p12s+p1s+p2s)/(24 m3s^2)),
s12,1/(12 m3s^2 (-1+t1)^4 t1) (-(-1+t1) (3 (p12s+p2s) t1 (5+t1)+p1s (-2-5 t1+t1^2))+6 t1 (p12s-p1s+p2s+2 p12s t1+2 p2s t1) Log[t1]),
s13,-1/(12 m3s^2 (-1+t2)^4) ((-1+t2) (-3 (p1s+p2s) (1+5 t2)+p12s (-1+5 t2+2 t2^2))+6 t2 (-p12s t2+p1s (2+t2)+p2s (2+t2)) Log[t2]),
s23,1/(12 m3s^2 (-1+t1)^4) ((-1+t1) (p2s-5 p2s t1-2 p2s t1^2+3 p12s (1+5 t1)+3 p1s (1+5 t1))-6 t1 (-p2s t1+p12s (2+t1)+p1s (2+t1)) Log[t1]),
us,(p12s (-t1 (t1+t1^2-2 t2) (-1+t2)^2 Log[t1]+(-1+t1) (-(-1+t2) (-t1^2 (-2+t2)+t1 (-3+t2) t2+t2^2)+(-1+t1)^2 t2^2 Log[t2])))/(2 m3s^2 (-1+t1)^3 (t1-t2)^2 (-1+t2)^2)-(p1s (t1 (-1+t2)^2 (t1^2+(-2+t1) t2) Log[t1]+(-1+t1) ((-1+t2) (t1^2-2 t1^2 t2+(-1+2 t1) t2^2)-(-1+t1) t2 (t1 (-2+t2)+t2^2) Log[t2])))/(2 m3s^2 (-1+t1)^2 (t1-t2)^3 (-1+t2)^2)+(p2s (t1^2 (-1+t2)^3 Log[t1]-(-1+t1) ((-1+t2) (2 t2^2+t1^2 (1+t2)-t1 t2 (3+t2))+t2 (-2 t1^2-t2 (1+t2)+t1 (2+t2+t2^2)) Log[t2])))/(2 m3s^2 (-1+t1)^2 (t1-t2)^2 (-1+t2)^3)
]];


CFant[LoopTools`cc1,m1s_,m2s_,m3s_]:=Block[{t1,t2,s123,s12,s13,s23,us},
t1=m1s/m3s;t2=m2s/m3s;
s123=SameQ[m1s,m2s,m3s];
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s];
s13=SameQ[m1s,m3s]&&UnsameQ[m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m1s,m3s];
us=UnsameQ[m1s,m2s,m3s];
Which[s123,1/(6m3s),
s12,(3-4t1+t1^2+2Log[t1])/(4m3s (-1+t1)^3),
s13,(-1+t2^2-2t2 Log[t2])/(2m3s (-1+t2)^3),
s23,(-1+4t1-3t1^2+2t1^2 Log[t1])/(4m3s (-1+t1)^3),
us,1/m3s (t1^2 (-1+t2)^2 Log[t1]-(-1+t1) t2 ((t1-t2) (-1+t2)+(t1 (-2+t2)+t2) Log[t2]))/(2(-1+t1) (t1-t2)^2 (-1+t2)^2)]];
CFDant[LoopTools`cc1,p1s_,p2s_,p12s_,m1s_,m2s_,m3s_]:=Block[{t1,t2,s123,s12,s13,s23,us},
t1=m1s/m3s;t2=m2s/m3s;
s123=SameQ[m1s,m2s,m3s];
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s];
s13=SameQ[m1s,m3s]&&UnsameQ[m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m1s,m3s];
us=UnsameQ[m1s,m2s,m3s];
Which[
s123,(p12s+2 (p1s+p2s))/(120 m3s^2),
s12,1/(36 m3s^2 (-1+t1)^5 t1) ((-1+t1) ((p12s+2 p2s) t1 (-17-8 t1+t1^2)+p1s (3+13 t1-5 t1^2+t1^3))+6 t1 (p12s-2 p1s+2 p2s+3 p12s t1+6 p2s t1) Log[t1]),
s13,-1/(36 m3s^2 (-1+t2)^5) ((-1+t2) (-6 (p1s+p2s) (1+10 t2+t2^2)+p12s (-1+8 t2+17 t2^2))+6 t2 (6 p1s (1+t2)+6 p2s (1+t2)-p12s t2 (3+t2)) Log[t2]),
s23,1/(36 m3s^2 (-1+t1)^5) ((-1+t1) (p12s (1-8 t1-17 t1^2)-2 p1s (-1+8 t1+17 t1^2)+p2s (1-5 t1+13 t1^2+3 t1^3))+6 t1^2 (-2 p2s t1+p12s (3+t1)+2 p1s (3+t1)) Log[t1]),
us,1/(6 m3s^2 (-1+t1)^3 (t1-t2)^3 (-1+t2)^3) p12s (t1^2 (-1+t2)^3 (t1+t1^2-3 t2+t1 t2) Log[t1]-(-1+t1) (2 (-1+t2) (-t2^3+t1 t2^2 (2+t2)-t1^2 t2 (2+t2^2)+t1^3 (1-t2+t2^2))+(-1+t1)^2 t2^2 (t1 (-3+t2)+t2 (1+t2)) Log[t2]))+1/(6 m3s^2 (-1+t1)^2 (t1-t2)^3 (-1+t2)^4) p2s (-2 t1^3 (-1+t2)^4 Log[t1]+(-1+t1) ((-1+t2) (t2^3 (5+t2)+t1^2 t2 (11+7 t2)+t1^3 (-2-5 t2+t2^2)-t1 t2^2 (14+3 t2+t2^2))+2 t2 (3 t1^3-t2^2 (1+2 t2)+t1 t2 (3+5 t2+t2^2)+t1^2 (-3-3 t2-4 t2^2+t2^3)) Log[t2]))+1/(6 m3s^2 (-1+t1)^2 (t1-t2)^4 (-1+t2)^3) p1s (2 t1^2 (-1+t2)^3 (t1^2-3 t2+2 t1 t2) Log[t1]-(-1+t1) ((-1+t2) (t2^3 (1+t2)+t1^2 t2 (3+7 t2-4 t2^2)+t1 t2^2 (-6+t2-t2^2)+t1^3 (2-9 t2+5 t2^2))+2 t2 (-t1 (-4+t2) t2^2-t2^3+t1^3 (3-3 t2+t2^2)+t1^2 (-3+3 t2-5 t2^2+2 t2^3)) Log[t2]))]];


CFant[LoopTools`cc2,m1s_,m2s_,m3s_]:=Block[{t1,t2,s123,s12,s13,s23,us},
t1=m1s/m3s;t2=m2s/m3s;
s123=SameQ[m1s,m2s,m3s];
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s];
s13=SameQ[m1s,m3s]&&UnsameQ[m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m1s,m3s];
us=UnsameQ[m1s,m2s,m3s];
Which[s123,1/(6m3s),
s12,(-1+t1^2-2t1 Log[t1])/(2m3s (-1+t1)^3),
s13,(-1+4t2-3t2^2+2t2^2 Log[t2])/(4m3s (-1+t2)^3),
s23,(-1+4t1-3t1^2+2t1^2 Log[t1])/(4m3s (-1+t1)^3),
us,(t1^2 (-1+t2)^2 Log[t1]-(-1+t1) (-(t1-t2) (-1+t2)+(-1+t1) t2^2 Log[t2]))/(2 m3s (-1+t1)^2 (t1-t2) (-1+t2)^2)]];
CFDant[LoopTools`cc2,p1s_,p2s_,p12s_,m1s_,m2s_,m3s_]:=Block[{t1,t2,s123,s12,s13,s23,us},
t1=m1s/m3s;t2=m2s/m3s;
s123=SameQ[m1s,m2s,m3s];
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s];
s13=SameQ[m1s,m3s]&&UnsameQ[m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m1s,m3s];
us=UnsameQ[m1s,m2s,m3s];
Which[
s123,(2 p12s+p1s+2 p2s)/(120 m3s^2),
s12,1/(36 m3s^2 (-1+t1)^5) ((-1+t1) (p1s (-17-8 t1+t1^2)+6 p12s (1+10 t1+t1^2)+6 p2s (1+10 t1+t1^2))+6 (p1s+3 p1s t1-6 (p12s+p2s) t1 (1+t1)) Log[t1]),
s13,1/(36 m3s^2 (-1+t2)^5) ((-1+t2) (-(p1s+2 p2s) (-1+8 t2+17 t2^2)+p12s (1-5 t2+13 t2^2+3 t2^3))+6 t2^2 (-2 p12s t2+p1s (3+t2)+2 p2s (3+t2)) Log[t2]),
s23,1/(36 m3s^2 (-1+t1)^5) ((-1+t1) (p1s+p2s-8 p1s t1-5 p2s t1-17 p1s t1^2+13 p2s t1^2+3 p2s t1^3-2 p12s (-1+8 t1+17 t1^2))+6 t1^2 (-2 p2s t1+2 p12s (3+t1)+p1s (3+t1)) Log[t1]),
us,1/(6 m3s^2 (-1+t1)^4 (t1-t2)^2 (-1+t2)^3) p12s (2 t1^2 (2 t1+t1^2-3 t2) (-1+t2)^3 Log[t1]-(-1+t1) ((-1+t2) (t2^2 (1+t2)+t1 t2 (-2+t2-5 t2^2)+t1^3 (5-9 t2+2 t2^2)+t1^2 (1-7 t2+14 t2^2-2 t2^3))+2 (-1+t1)^3 t2^3 Log[t2]))+1/(6 m3s^2 (-1+t1)^3 (t1-t2)^3 (-1+t2)^3) p1s (t1^2 (-1+t2)^3 (t1+t1^2-3 t2+t1 t2) Log[t1]-(-1+t1) (2 (-1+t2) (-t2^3+t1 t2^2 (2+t2)-t1^2 t2 (2+t2^2)+t1^3 (1-t2+t2^2))+(-1+t1)^2 t2^2 (t1 (-3+t2)+t2 (1+t2)) Log[t2]))+1/(6 m3s^2 (-1+t1)^3 (t1-t2)^2 (-1+t2)^4) p2s (-2 t1^3 (-1+t2)^4 Log[t1]+(-1+t1) (-(-1+t2) (t2^2 (1+5 t2)+t1^3 (1-5 t2-2 t2^2)-t1 t2 (2+7 t2+9 t2^2)+t1^2 (1+t2+14 t2^2+2 t2^3))-2 (-1+t1)^2 t2^2 (3 t1-t2 (2+t2)) Log[t2]))]];



CFant[LoopTools`cc00,m1s_,m2s_,m3s_]:=Block[{t1,t2,s123,s12,s13,s23,us},
t1=m1s/m3s;t2=m2s/m3s;
s123=SameQ[m1s,m2s,m3s];
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s];
s13=SameQ[m1s,m3s]&&UnsameQ[m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m1s,m3s];
us=UnsameQ[m1s,m2s,m3s];
Which[s123,-(Log[m3s]/4),
s12,(3-4 t1+t1^2-2 (-1+t1)^2 Log[m3s]-2 (-2+t1) t1 Log[t1])/(8 (-1+t1)^2),
s13,(1-4 t2+3 t2^2-2 (-1+t2)^2 Log[m3s]-2 t2^2 Log[t2])/(8 (-1+t2)^2),
s23,(1-4 t1+3 t1^2-2 (-1+t1)^2 Log[m3s]-2 t1^2 Log[t1])/(8 (-1+t1)^2),
us,-(1/4)Log[m3s]+(-2 t1^2 (-1+t2) Log[t1]+(-1+t1) (3 (t1-t2) (-1+t2)+2 t2^2 Log[t2]))/(8 (-1+t1) (t1-t2) (-1+t2))]];
CFDant[LoopTools`cc00,p1s_,p2s_,p12s_,m1s_,m2s_,m3s_]:=Block[{t1,t2,s123,s12,s13,s23,us},
t1=m1s/m3s;t2=m2s/m3s;
s123=SameQ[m1s,m2s,m3s];
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s];
s13=SameQ[m1s,m3s]&&UnsameQ[m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m1s,m3s];
us=UnsameQ[m1s,m2s,m3s];
Which[s123,(p12s+p1s+p2s)/(48 m3s),
s12,1/(72 m3s (-1+t1)^4) ((-1+t1) (3 p12s (-2-5 t1+t1^2)+3 p2s (-2-5 t1+t1^2)+p1s (11-7 t1+2 t1^2))-6 (p1s-3 (p12s+p2s) t1) Log[t1]),
s13,(-(-1+t2) (-3 (p1s+p2s) (-1+5 t2+2 t2^2)+p12s (2-7 t2+11 t2^2))+6 t2^2 (-3 p1s-3 p2s+p12s t2) Log[t2])/(72 m3s (-1+t2)^4),
s23,1/(72 m3s (-1+t1)^4) ((-1+t1) (p2s (-2+7 t1-11 t1^2)+3 p12s (-1+5 t1+2 t1^2)+3 p1s (-1+5 t1+2 t1^2))+6 t1^2 (-3 p12s-3 p1s+p2s t1) Log[t1]),
us,(p12s (-t1^2 (-1+t2)^2 (-3 t2+t1 (2+t2)) Log[t1]+(-1+t1) (-(-1+t2) (-t1^3 (-1+t2)-2 t1 t2+t2^2+t1^2 (1-t2+t2^2))+(-1+t1)^2 t2^3 Log[t2])))/(12 m3s (-1+t1)^3 (t1-t2)^2 (-1+t2)^2)+(p2s (t1^3 (-1+t2)^3 Log[t1]-(-1+t1) ((-1+t2) (t2^2 (1+t2)+t1^2 (1+t2^2)-t1 t2 (2+t2+t2^2))+(-1+t1) t2^2 (t1 (-3+t2)+2 t2) Log[t2])))/(12 m3s (-1+t1)^2 (t1-t2)^2 (-1+t2)^3)+(p1s (-t1^2 (-1+t2)^2 (t1-3 t2+2 t1 t2) Log[t1]+(-1+t1) ((-1+t2) (t1^3 (-1+t2)+t1^2 t2+t2^3-t1 t2^2 (1+t2))+(-1+t1) t2^2 (t2+t1 (-3+2 t2)) Log[t2])))/(12 m3s (-1+t1)^2 (t1-t2)^3 (-1+t2)^2)
]];


CFant[LoopTools`cc11,m1s_,m2s_,m3s_]:=Block[{t1,t2,s123,s12,s13,s23,us},
t1=m1s/m3s;t2=m2s/m3s;
s123=SameQ[m1s,m2s,m3s];
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s];
s13=SameQ[m1s,m3s]&&UnsameQ[m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m1s,m3s];
us=UnsameQ[m1s,m2s,m3s];
Which[s123,-(1/(12 m3s)),
s12,(11-18 t1+9 t1^2-2 t1^3+6 Log[t1])/(18 m3s (-1+t1)^4),
s13,-((2+3 t2-6 t2^2+t2^3+6 t2 Log[t2])/(6 m3s (-1+t2)^4)),
s23,(-2+9 t1-18 t1^2+11 t1^3-6 t1^3 Log[t1])/(18 m3s (-1+t1)^4),
us,
1/m3s 1/(6(-1+t1) (t1-t2)^3 (-1+t2)^3) (-2 t1^3 (-1+t2)^3 Log[t1]+(-1+t1) t2 ((-1+t2) (-4 t1 (-2+t2) t2+(-3+t2) t2^2+t1^2 (-5+3 t2))+2 (t1 (-3+t2) t2+t2^2+t1^2 (3-3 t2+t2^2)) Log[t2]))]
];
CFDant[LoopTools`cc11,p1s_,p2s_,p12s_,m1s_,m2s_,m3s_]:=-(1/(48 (m1s-m2s)^4 (m1s-m3s)^2 (m2s-m3s)^4))((m1s-m2s) (-3 m2s^3 m3s^2 (m2s^2-4 m2s m3s-13 m3s^2)+3 m1s m2s^2 m3s (2 m2s^3-7 m2s^2 m3s-22 m2s m3s^2-37 m3s^3)+4 m1s^4 (m2s^3-5 m2s^2 m3s+13 m2s m3s^2+3 m3s^3)+m1s^3 (5 m2s^4-10 m2s^3 m3s-41 m2s^2 m3s^2-130 m2s m3s^3-16 m3s^4)+m1s^2 m2s (-3 m2s^4+49 m2s^2 m3s^2+148 m2s m3s^3+94 m3s^4)) p12s+(m1s-m3s) ((m2s^3 m3s (-3 m2s^3+22 m2s^2 m3s+31 m2s m3s^2-2 m3s^3)+4 m1s^4 (m2s^3-5 m2s^2 m3s+13 m2s m3s^2+3 m3s^3)+m1s^2 m2s (-20 m2s^4+59 m2s^3 m3s+42 m2s^2 m3s^2+233 m2s m3s^3-26 m3s^4)-m1s^3 (11 m2s^4-62 m2s^3 m3s+187 m2s^2 m3s^2+50 m2s m3s^3+6 m3s^4)+m1s m2s^2 (3 m2s^4-2 m2s^3 m3s-73 m2s^2 m3s^2-130 m2s m3s^3+10 m3s^4)) p1s-(m1s-m2s) (m2s^3 m3s (m2s^2-8 m2s m3s-17 m3s^2)+m1s^2 m2s (5 m2s^3-12 m2s^2 m3s-21 m2s m3s^2-44 m3s^3)+2 m1s^3 (m2s^3-5 m2s^2 m3s+13 m2s m3s^2+3 m3s^3)+m1s m2s^2 (-m2s^3+3 m2s^2 m3s+21 m2s m3s^2+49 m3s^3)) p2s))+1/(24 (m1s-m2s)^5 (m2s-m3s)^5) m2s ((m1s-m2s) (24 m1s^3 m3s^3+m2s^3 m3s (2 m2s^2-19 m2s m3s-7 m3s^2)+2 m1s m2s^2 (m2s^3-9 m2s^2 m3s+30 m2s m3s^2+14 m3s^3)-m1s^2 m2s (5 m2s^3-25 m2s^2 m3s+50 m2s m3s^2+42 m3s^3)) p12s+3 ((8 m1s^4 m3s^3+m2s^4 m3s^2 (7 m2s+m3s)+m1s m2s^3 m3s (6 m2s^2-33 m2s m3s-5 m3s^2)+m1s^2 m2s^2 (5 m2s^3-27 m2s^2 m3s+60 m2s m3s^2+10 m3s^3)-m1s^3 (m2s^4-5 m2s^3 m3s+10 m2s^2 m3s^2+30 m2s m3s^3-4 m3s^4)) p1s-(m1s-m2s) (4 m1s^3 m3s^3-m2s^3 m3s^2 (3 m2s+m3s)-2 m1s m2s^2 m3s (m2s^2-5 m2s m3s-2 m3s^2)-m1s^2 m2s (m2s^3-5 m2s^2 m3s+10 m2s m3s^2+6 m3s^3)) p2s)) Log[m2s/m1s]+(m3s^3 ((24 m1s^2 m2s+m1s (4 m2s^2-47 m2s m3s-5 m3s^2)+m3s (-2 m2s^2+19 m2s m3s+7 m3s^2)) p12s+3 (m1s-m3s) ((8 m1s m2s-m3s (7 m2s+m3s)) p1s+(-4 m1s m2s+m3s (3 m2s+m3s)) p2s)) Log[m3s/m1s])/(24 (m1s-m3s)^3 (-m2s+m3s)^5);


CFant[LoopTools`cc12,m1s_,m2s_,m3s_]:=Block[{t1,t2,s123,s12,s13,s23,us},
t1=m1s/m3s;t2=m2s/m3s;
s123=SameQ[m1s,m2s,m3s];
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s];
s13=SameQ[m1s,m3s]&&UnsameQ[m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m1s,m3s];
us=UnsameQ[m1s,m2s,m3s];
Which[s123,-(1/(24 m3s)),
s12,-((2+3 t1-6 t1^2+t1^3+6 t1 Log[t1])/(12 m3s (-1+t1)^4)),
s13,(-1+6 t2-3 t2^2-2 t2^3+6 t2^2 Log[t2])/(12 m3s (-1+t2)^4),
s23,(-2+9 t1-18 t1^2+11 t1^3-6 t1^3 Log[t1])/(36 m3s (-1+t1)^4),
us,
1/m3s (-t1^3 (-1+t2)^3 Log[t1]+(-1+t1) ((-1+t2) (t2^2 (1+t2)+t1^2 (1+t2^2)-t1 t2 (2+t2+t2^2))+(-1+t1) t2^2 (t1 (-3+t2)+2 t2) Log[t2]))/(6 (-1+t1)^2 (t1-t2)^2 (-1+t2)^3)]];
CFDant[LoopTools`cc12,p1s_,p2s_,p12s_,m1s_,m2s_,m3s_]:=1/(48 (m1s-m2s)^3 (m1s-m3s)^3 (m2s-m3s)^4) ((m2s^3 m3s^3 (m2s^2-44 m2s m3s-5 m3s^2)-4 m1s^5 (m2s^3-7 m2s^2 m3s-7 m2s m3s^2+m3s^3)+m1s m2s^2 m3s^2 (-m2s^3+113 m2s^2 m3s+113 m2s m3s^2+15 m3s^3)+m1s^4 (4 m2s^4-46 m2s^3 m3s-149 m2s^2 m3s^2-54 m2s m3s^3+5 m3s^4)+m1s^3 m3s (24 m2s^4+189 m2s^3 m3s+251 m2s^2 m3s^2+11 m2s m3s^3+5 m3s^4)-m1s^2 m2s m3s (6 m2s^4+67 m2s^3 m3s+307 m2s^2 m3s^2+85 m2s m3s^3+15 m3s^4)) p12s-(m1s-m3s) ((m2s^3 m3s^2 (-5 m2s^2-44 m2s m3s+m3s^2)+m1s m2s^2 m3s (10 m2s^3+69 m2s^2 m3s+114 m2s m3s^2-m3s^3)+4 m1s^4 (m2s^3-7 m2s^2 m3s-7 m2s m3s^2+m3s^3)+m1s^3 m2s (-5 m2s^3+58 m2s^2 m3s+121 m2s m3s^2+18 m3s^3)-m1s^2 m2s (5 m2s^4+16 m2s^3 m3s+193 m2s^2 m3s^2+68 m2s m3s^3+6 m3s^4)) p1s-2 (-m2s^3 m3s^2 (m2s^2+10 m2s m3s+m3s^2)+m1s^4 (m2s^3-7 m2s^2 m3s-7 m2s m3s^2+m3s^3)+m1s^3 m3s (10 m2s^3+35 m2s^2 m3s+2 m2s m3s^2+m3s^3)+m1s m2s^2 m3s (2 m2s^3+17 m2s^2 m3s+26 m2s m3s^2+3 m3s^3)-m1s^2 m2s (m2s^4+5 m2s^3 m3s+44 m2s^2 m3s^2+19 m2s m3s^3+3 m3s^4)) p2s))-1/(24 (m1s-m2s)^4 (m2s-m3s)^5) m2s^2 ((m1s-m2s) (24 m1s^2 m3s^2+3 m1s m2s (m2s^2-5 m2s m3s-12 m3s^2)+m2s^2 (-m2s^2+11 m2s m3s+14 m3s^2)) p12s+(24 m1s^3 m3s^2+m2s^2 m3s (-14 m2s^2-11 m2s m3s+m3s^2)+3 m1s^2 (m2s^3-5 m2s^2 m3s-22 m2s m3s^2+2 m3s^3)+m1s (-6 m2s^4+38 m2s^3 m3s+44 m2s^2 m3s^2-4 m2s m3s^3)) p1s-2 (m1s-m2s) (6 m1s^2 m3s^2+3 m2s^2 m3s (m2s+m3s)+m1s m2s (m2s^2-5 m2s m3s-8 m3s^2)) p2s) Log[m2s/m1s]-1/(24 (m1s-m3s)^4 (-m2s+m3s)^5) m3s^2 ((24 m1s^3 m2s^2+m2s m3s^2 (m2s^2-11 m2s m3s-14 m3s^2)+3 m1s^2 (2 m2s^3-22 m2s^2 m3s-5 m2s m3s^2+m3s^3)+m1s (-4 m2s^3 m3s+44 m2s^2 m3s^2+38 m2s m3s^3-6 m3s^4)) p12s+(m1s-m3s) ((24 m1s^2 m2s^2+m3s^2 (14 m2s^2+11 m2s m3s-m3s^2)+3 m1s m3s (-12 m2s^2-5 m2s m3s+m3s^2)) p1s-2 (6 m1s^2 m2s^2+3 m2s m3s^2 (m2s+m3s)+m1s m3s (-8 m2s^2-5 m2s m3s+m3s^2)) p2s)) Log[m3s/m1s];


CFant[LoopTools`cc22,m1s_,m2s_,m3s_]:=Block[{t1,t2,s123,s12,s13,s23,us},
t1=m1s/m3s;t2=m2s/m3s;
s123=SameQ[m1s,m2s,m3s];
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s];
s13=SameQ[m1s,m3s]&&UnsameQ[m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m1s,m3s];
us=UnsameQ[m1s,m2s,m3s];
Which[s123,-(1/(12 m3s)),
s12,(-1+6 t1-3 t1^2-2 t1^3+6 t1^2 Log[t1])/(6 m3s (-1+t1)^4),
s13,(-2+9 t2-18 t2^2+11 t2^3-6 t2^3 Log[t2])/(18 m3s (-1+t2)^4),
s23,(-2+9 t1-18 t1^2+11 t1^3-6 t1^3 Log[t1])/(18 m3s (-1+t1)^4),
us,
1/m3s (-2 t1^3 (-1+t2)^3 Log[t1]+(-1+t1) ((-1+t2) (t1^2 (3-5 t2)+t2-3 t2^2+t1 (-1+5 t2^2))+2 (-1+t1)^2 t2^3 Log[t2]))/(6 (-1+t1)^3 (t1-t2) (-1+t2)^3)]];
CFDant[LoopTools`cc22,p1s_,p2s_,p12s_,m1s_,m2s_,m3s_]:=-(1/(48 (m1s-m2s)^2 (m1s-m3s)^4 (m2s-m3s)^4))((4 m1s^5 (3 m2s^3+13 m2s^2 m3s-5 m2s m3s^2+m3s^3)+m2s^2 m3s^3 (2 m2s^3-31 m2s^2 m3s-22 m2s m3s^2+3 m3s^3)-2 m1s m2s m3s^2 (5 m2s^4-64 m2s^3 m3s-52 m2s^2 m3s^2-12 m2s m3s^3+3 m3s^4)-m1s^4 (18 m2s^4+102 m2s^3 m3s+167 m2s^2 m3s^2-58 m2s m3s^3+11 m3s^4)+m1s^3 (6 m2s^5+24 m2s^4 m3s+420 m2s^3 m3s^2-20 m2s^2 m3s^3+70 m2s m3s^4-20 m3s^5)+m1s^2 m3s (26 m2s^5-223 m2s^4 m3s-172 m2s^3 m3s^2-132 m2s^2 m3s^3+18 m2s m3s^4+3 m3s^5)) p12s+(m1s-m3s) ((3 m2s^2 m3s^3 (13 m2s^2+4 m2s m3s-m3s^2)-3 m1s m2s m3s^2 (37 m2s^3+22 m2s^2 m3s+7 m2s m3s^2-2 m3s^3)+4 m1s^4 (3 m2s^3+13 m2s^2 m3s-5 m2s m3s^2+m3s^3)-m1s^3 (16 m2s^4+130 m2s^3 m3s+41 m2s^2 m3s^2+10 m2s m3s^3-5 m3s^4)+m1s^2 m3s (94 m2s^4+148 m2s^3 m3s+49 m2s^2 m3s^2-3 m3s^4)) p1s+(m2s^2 m3s^3 (-17 m2s^2-8 m2s m3s+m3s^2)+m1s m2s m3s^2 (49 m2s^3+38 m2s^2 m3s+11 m2s m3s^2-2 m3s^3)-2 m1s^4 (3 m2s^3+13 m2s^2 m3s-5 m2s m3s^2+m3s^3)+m1s^3 (6 m2s^4+70 m2s^3 m3s+11 m2s^2 m3s^2+14 m2s m3s^3-5 m3s^4)+m1s^2 m3s (-44 m2s^4-70 m2s^3 m3s-33 m2s^2 m3s^2+2 m2s m3s^3+m3s^4)) p2s))+(m2s^3 (3 (m1s-m2s) (-m2s^2+8 m1s m3s-7 m2s m3s) p12s+(-5 m1s m2s^2+7 m2s^3+24 m1s^2 m3s-47 m1s m2s m3s+19 m2s^2 m3s+4 m1s m3s^2-2 m2s m3s^2) p1s-3 (m1s-m2s) (-m2s^2+4 m1s m3s-3 m2s m3s) p2s) Log[m2s/m1s])/(24 (m1s-m2s)^3 (m2s-m3s)^5)+1/(24 (m1s-m3s)^5 (-m2s+m3s)^5) m3s (3 (8 m1s^4 m2s^3+m2s^2 m3s^4 (m2s+7 m3s)+m1s m2s m3s^3 (-5 m2s^2-33 m2s m3s+6 m3s^2)+m1s^2 m3s^2 (10 m2s^3+60 m2s^2 m3s-27 m2s m3s^2+5 m3s^3)+m1s^3 (4 m2s^4-30 m2s^3 m3s-10 m2s^2 m3s^2+5 m2s m3s^3-m3s^4)) p12s+(m1s-m3s) ((24 m1s^3 m2s^3+m2s m3s^3 (-7 m2s^2-19 m2s m3s+2 m3s^2)+2 m1s m3s^2 (14 m2s^3+30 m2s^2 m3s-9 m2s m3s^2+m3s^3)-m1s^2 m3s (42 m2s^3+50 m2s^2 m3s-25 m2s m3s^2+5 m3s^3)) p1s+3 (-4 m1s^3 m2s^3+m2s^2 m3s^3 (m2s+3 m3s)-2 m1s m2s m3s^2 (2 m2s^2+5 m2s m3s-m3s^2)+m1s^2 m3s (6 m2s^3+10 m2s^2 m3s-5 m2s m3s^2+m3s^3)) p2s)) Log[m3s/m1s];


CFDant[id_,p1s_,p2s_,p12s_,m1s_,m2s_,m3s_]:=Block[{},
Print["Derivative for C0i[",id,"] not implemented"];
C0i[id,p1s,p2s,p12s,m1s,m2s,m3s]-CFant[id,m1s,m2s,m3s]];


(* ::Subsection::Closed:: *)
(*D*)


DFant[LoopTools`dd0,m1s_,m2s_,m3s_,m4s_]:=Block[{t1,t2,t3,s12,s13,s14,s23,s24,s34,s12s34,s13s24,s14s23,s123,s124,s134,s234,s1234,us},
{t1,t2,t3}={m1s/m4s,m2s/m4s,m3s/m4s};
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s,m4s];
s13=SameQ[m1s,m3s]&&UnsameQ[m1s,m2s,m4s];
s14=SameQ[m1s,m4s]&&UnsameQ[m1s,m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m2s,m1s,m4s];
s24=SameQ[m2s,m4s]&&UnsameQ[m2s,m1s,m3s];
s34=SameQ[m3s,m4s]&&UnsameQ[m3s,m1s,m2s];
s12s34=SameQ[m1s,m2s]&&SameQ[m3s,m4s]&&UnsameQ[m1s,m3s];
s13s24=SameQ[m1s,m3s]&&SameQ[m2s,m4s]&&UnsameQ[m1s,m2s];
s14s23=SameQ[m1s,m4s]&&SameQ[m2s,m3s]&&UnsameQ[m1s,m2s];
s123=SameQ[m1s,m2s,m3s]&&UnsameQ[m1s,m4s];
s124=SameQ[m1s,m2s,m4s]&&UnsameQ[m1s,m3s];
s134=SameQ[m1s,m3s,m4s]&&UnsameQ[m1s,m2s];
s234=SameQ[m2s,m3s,m4s]&&UnsameQ[m2s,m1s];
s1234=SameQ[m1s,m2s,m3s,m4s];
us=UnsameQ[m1s,m2s,m3s,m4s];

Which[s12,((t1^2-t3) (-1+t3) Log[t1]-(-1+t1) ((t1-t3) (-1+t3)+(-1+t1) t3 Log[t3]))/(m4s^2 (-1+t1)^2 (t1-t3)^2 (-1+t3)),s13,((t1^2-t2) (-1+t2) Log[t1]-(-1+t1) ((t1-t2) (-1+t2)+(-1+t1) t2 Log[t2]))/(m4s^2 (-1+t1)^2 (t1-t2)^2 (-1+t2)),
s14,(-t2 (-1+t3)^2 Log[t2]+(-1+t2) (-(t2-t3) (-1+t3)+(-1+t2) t3 Log[t3]))/(m4s^2 (-1+t2)^2 (t2-t3) (-1+t3)^2),
s23,(-t1 (-1+t2)^2 Log[t1]+(-1+t1) ((t1-t2) (-1+t2)+(-t1+t2^2) Log[t2]))/(m4s^2 (-1+t1) (t1-t2)^2 (-1+t2)^2),
s24,(-t1 (-1+t3)^2 Log[t1]+(-1+t1) (-(t1-t3) (-1+t3)+(-1+t1) t3 Log[t3]))/(m4s^2 (-1+t1)^2 (t1-t3) (-1+t3)^2),
s34,(-t1 (-1+t2)^2 Log[t1]+(-1+t1) (-(t1-t2) (-1+t2)+(-1+t1) t2 Log[t2]))/(m4s^2 (-1+t1)^2 (t1-t2) (-1+t2)^2),
s12s34,(2-2 t1+(1+t1) Log[t1])/(m4s^2 (-1+t1)^3),
s13s24,(2-2 t1+(1+t1) Log[t1])/(m4s^2 (-1+t1)^3),
s14s23,(2-2 t2+(1+t2) Log[t2])/(m4s^2 (-1+t2)^3),
s123,(-1+t1^2-2 t1 Log[t1])/(2 m4s^2 (-1+t1)^3 t1),
s124,(-1+t3^2-2 t3 Log[t3])/(2 m4s^2 (-1+t3)^3),
s134,(-1+t2^2-2 t2 Log[t2])/(2 m4s^2 (-1+t2)^3),
s234,(-1+t1^2-2 t1 Log[t1])/(2 m4s^2 (-1+t1)^3),
s1234,1/(6 m4s^2),
us,(-t1 (-1+t2) (t2-t3) (-1+t3) Log[t1]+(-1+t1) (t2 (t1-t3) (-1+t3) Log[t2]-(t1-t2) (-1+t2) t3 Log[t3]))/(m4s^2 (-1+t1) (t1-t2) (-1+t2) (t1-t3) (t2-t3) (-1+t3))]];


DFant[LoopTools`dd1,m1s_,m2s_,m3s_,m4s_]:=Block[{t1,t2,t3,s12,s13,s14,s23,s24,s34,s12s34,s13s24,s14s23,s123,s124,s134,s234,s1234,us},
{t1,t2,t3}={m1s/m4s,m2s/m4s,m3s/m4s};
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s,m4s];
s13=SameQ[m1s,m3s]&&UnsameQ[m1s,m2s,m4s];
s14=SameQ[m1s,m4s]&&UnsameQ[m1s,m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m2s,m1s,m4s];
s24=SameQ[m2s,m4s]&&UnsameQ[m2s,m1s,m3s];
s34=SameQ[m3s,m4s]&&UnsameQ[m3s,m1s,m2s];
s12s34=SameQ[m1s,m2s]&&SameQ[m3s,m4s]&&UnsameQ[m1s,m3s];
s13s24=SameQ[m1s,m3s]&&SameQ[m2s,m4s]&&UnsameQ[m1s,m2s];
s14s23=SameQ[m1s,m4s]&&SameQ[m2s,m3s]&&UnsameQ[m1s,m2s];
s123=SameQ[m1s,m2s,m3s]&&UnsameQ[m1s,m4s];
s124=SameQ[m1s,m2s,m4s]&&UnsameQ[m1s,m3s];
s134=SameQ[m1s,m3s,m4s]&&UnsameQ[m1s,m2s];
s234=SameQ[m2s,m3s,m4s]&&UnsameQ[m2s,m1s];
s1234=SameQ[m1s,m2s,m3s,m4s];
us=UnsameQ[m1s,m2s,m3s,m4s];
Which[s12,(2 (-1+t3) (-3 t1^2 t3+t3^2+t1^3 (1+t3)) Log[t1]-(-1+t1) ((-1+t3) (t1^2+t1^3+3 t3^2-t1 t3 (4+t3))+2 (-1+t1)^2 t3^2 Log[t3]))/(4 m4s^2 (-1+t1)^3 (t1-t3)^3 (-1+t3)),
s13,-(t1 (-1+t2)^2 (t1^2-2 t2+t1 t2) Log[t1]-(-1+t1) ((-1+t2) (t2^2-2 t1 t2^2+t1^2 (-1+2 t2))+t2 (t1^2 (-2+t2)-t2^2+t1 (2-t2+t2^2)) Log[t2]))/(2 m4s^2 (-1+t1)^2 (t1-t2)^3 (-1+t2)^2),s14,(-t2 (t2+t2^2-2 t3) (-1+t3)^2 Log[t2]+(-1+t2) (-(-1+t3) (-t2^2 (-2+t3)+t2 (-3+t3) t3+t3^2)+(-1+t2)^2 t3^2 Log[t3]))/(2 m4s^2 (-1+t2)^3 (t2-t3)^2 (-1+t3)^2),
s23,(2 t1^2 (-1+t2)^3 Log[t1]-(-1+t1) (-(-1+t2) (-t1^2 (-3+t2)-4 t1 t2+t2^2 (1+t2))+2 (t1^2+t1 (-3+t2) t2^2+t2^3) Log[t2]))/(4 m4s^2 (-1+t1) (t1-t2)^3 (-1+t2)^3),
s24,(2 t1^2 (-1+t3)^3 Log[t1]-(-1+t1) ((-1+t3) (t1+t1^2 (1-3 t3)+3 t1 t3^2-t3 (1+t3))+2 (-1+t1)^2 t3^2 Log[t3]))/(4 m4s^2 (-1+t1)^3 (t1-t3) (-1+t3)^3),
s34,(t1^2 (-1+t2)^3 Log[t1]-(-1+t1) ((-1+t2) (2 t2^2+t1^2 (1+t2)-t1 t2 (3+t2))+t2 (-2 t1^2-t2 (1+t2)+t1 (2+t2+t2^2)) Log[t2]))/(2 m4s^2 (-1+t1)^2 (t1-t2)^2 (-1+t2)^3),
s12s34,(5-4 t1-t1^2+(2+4 t1) Log[t1])/(4 m4s^2 (-1+t1)^4),
s13s24,-((1+4 t1-5 t1^2+2 t1 (2+t1) Log[t1])/(4 m4s^2 (-1+t1)^4)),
s14s23,(2 t1^2 (-1+t2)^3 Log[t1]-(-1+t1) (-(-1+t2) (-t1^2 (-3+t2)-4 t1 t2+t2^2 (1+t2))+2 (t1^2+t1 (-3+t2) t2^2+t2^3) Log[t2]))/(4 m4s^2 (-1+t1) (t1-t2)^3 (-1+t2)^3),
s123,-((3+2/t1-6 t1+t1^2+6 Log[t1])/(12 m4s^2 (-1+t1)^4)),
s124,(-1+6 t3-3 t3^2-2 t3^3+6 t3^2 Log[t3])/(12 m4s^2 (-1+t3)^4),
s134,-((1+4 t2-5 t2^2+2 t2 (2+t2) Log[t2])/(4 m4s^2 (-1+t2)^4)),
s234,(-1+6 t1-3 t1^2-2 t1^3+6 t1^2 Log[t1])/(12 m4s^2 (-1+t1)^4),s1234,-1/(24 m4s^2),
us,(t1^2 (-1+t2)^2 (t2-t3)^2 (-1+t3) Log[t1]-(-1+t1) (t2 (-1+t3) (t2 t3 (-t2^2+t3)-t1^2 (t2-2 t3+t2 t3)+t1 (t2^3-2 t3^2+t2 t3^2)) Log[t2]+(t1-t2) (-1+t2) (t2 (t1-t3) (t2-t3) (-1+t3)+(t1-t2) (-1+t2) t3^2 Log[t3])))/(2 m4s^2 (-1+t1) (t1-t2)^2 (-1+t2)^2 (t1-t3) (t2-t3)^2 (-1+t3))]];


DFant[LoopTools`dd2,m1s_,m2s_,m3s_,m4s_]:=Block[{t1,t2,t3,s12,s13,s14,s23,s24,s34,s12s34,s13s24,s14s23,s123,s124,s134,s234,s1234,us},
{t1,t2,t3}={m1s/m4s,m2s/m4s,m3s/m4s};
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s,m4s];
s13=SameQ[m1s,m3s]&&UnsameQ[m1s,m2s,m4s];
s14=SameQ[m1s,m4s]&&UnsameQ[m1s,m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m2s,m1s,m4s];
s24=SameQ[m2s,m4s]&&UnsameQ[m2s,m1s,m3s];
s34=SameQ[m3s,m4s]&&UnsameQ[m3s,m1s,m2s];
s12s34=SameQ[m1s,m2s]&&SameQ[m3s,m4s]&&UnsameQ[m1s,m3s];
s13s24=SameQ[m1s,m3s]&&SameQ[m2s,m4s]&&UnsameQ[m1s,m2s];
s14s23=SameQ[m1s,m4s]&&SameQ[m2s,m3s]&&UnsameQ[m1s,m2s];
s123=SameQ[m1s,m2s,m3s]&&UnsameQ[m1s,m4s];
s124=SameQ[m1s,m2s,m4s]&&UnsameQ[m1s,m3s];
s134=SameQ[m1s,m3s,m4s]&&UnsameQ[m1s,m2s];
s234=SameQ[m2s,m3s,m4s]&&UnsameQ[m2s,m1s];
s1234=SameQ[m1s,m2s,m3s,m4s];
us=UnsameQ[m1s,m2s,m3s,m4s];
Which[s12,-(t1 (-1+t3)^2 (t1^2-2 t3+t1 t3) Log[t1]-(-1+t1) ((-1+t3) (t3^2-2 t1 t3^2+t1^2 (-1+2 t3))+t3 (t1^2 (-2+t3)-t3^2+t1 (2-t3+t3^2)) Log[t3]))/(2 m4s^2 (-1+t1)^2 (t1-t3)^3 (-1+t3)^2),
s13,(2 (-1+t2) (-3 t1^2 t2+t2^2+t1^3 (1+t2)) Log[t1]-(-1+t1) ((-1+t2) (t1^2+t1^3+3 t2^2-t1 t2 (4+t2))+2 (-1+t1)^2 t2^2 Log[t2]))/(4 m4s^2 (-1+t1)^3 (t1-t2)^3 (-1+t2)),
s14,(t2^2 (-1+t3)^3 Log[t2]-(-1+t2) ((-1+t3) (2 t3^2+t2^2 (1+t3)-t2 t3 (3+t3))+t3 (-2 t2^2-t3 (1+t3)+t2 (2+t3+t3^2)) Log[t3]))/(2 m4s^2 (-1+t2)^2 (t2-t3)^2 (-1+t3)^3),
s23,(2 t1^2 (-1+t2)^3 Log[t1]-(-1+t1) (-(-1+t2) (-t1^2 (-3+t2)-4 t1 t2+t2^2 (1+t2))+2 (t1^2+t1 (-3+t2) t2^2+t2^3) Log[t2]))/(4 m4s^2 (-1+t1) (t1-t2)^3 (-1+t2)^3),
s24,(t1^2 (-1+t3)^3 Log[t1]-(-1+t1) ((-1+t3) (2 t3^2+t1^2 (1+t3)-t1 t3 (3+t3))+t3 (-2 t1^2-t3 (1+t3)+t1 (2+t3+t3^2)) Log[t3]))/(2 m4s^2 (-1+t1)^2 (t1-t3)^2 (-1+t3)^3),
s34,(t1^2 (-1+t2)^3 Log[t1]+(-1+t1) (-(1/2) (t1-t2) (-1+t2) (1+t1+t2-3 t1 t2)-(-1+t1)^2 t2^2 Log[t2]))/(2 m4s^2 (-1+t1)^3 (t1-t2) (-1+t2)^3),
s12s34,-((1+4 t1-5 t1^2+2 t1 (2+t1) Log[t1])/(4 m4s^2 (-1+t1)^4)),
s13s24,(5-4 t1-t1^2+(2+4 t1) Log[t1])/(4 m4s^2 (-1+t1)^4),
s14s23,(5-4 t2-t2^2+(2+4 t2) Log[t2])/(4 m4s^2 (-1+t2)^4),
s123,-((3+2/t1-6 t1+t1^2+6 Log[t1])/(12 m4s^2 (-1+t1)^4)),
s124,-((1+4 t3-5 t3^2+2 t3 (2+t3) Log[t3])/(4 m4s^2 (-1+t3)^4)),
s134,(-1+6 t2-3 t2^2-2 t2^3+6 t2^2 Log[t2])/(12 m4s^2 (-1+t2)^4),
s234,(-1+6 t1-3 t1^2-2 t1^3+6 t1^2 Log[t1])/(12 m4s^2 (-1+t1)^4),
s1234,-(1/(24 m4s^2)),
us,-(-t1^2 (-1+t2) (t2-t3)^2 (-1+t3)^2 Log[t1]+(-1+t1) (t2^2 (t1-t3)^2 (-1+t3)^2 Log[t2]-(t1-t2) (-1+t2) t3 ((t1-t3) (t2-t3) (-1+t3)+(t1 (t2 (-2+t3)+t3)+t3 (t2-t3^2)) Log[t3])))/(2 m4s^2 (-1+t1) (t1-t2) (-1+t2) (t1-t3)^2 (t2-t3)^2 (-1+t3)^2)]];


DFant[LoopTools`dd3,m1s_,m2s_,m3s_,m4s_]:=Block[{t1,t2,t3,s12,s13,s14,s23,s24,s34,s12s34,s13s24,s14s23,s123,s124,s134,s234,s1234,us},
{t1,t2,t3}={m1s/m4s,m2s/m4s,m3s/m4s};
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s,m4s];
s13=SameQ[m1s,m3s]&&UnsameQ[m1s,m2s,m4s];
s14=SameQ[m1s,m4s]&&UnsameQ[m1s,m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m2s,m1s,m4s];
s24=SameQ[m2s,m4s]&&UnsameQ[m2s,m1s,m3s];
s34=SameQ[m3s,m4s]&&UnsameQ[m3s,m1s,m2s];
s12s34=SameQ[m1s,m2s]&&SameQ[m3s,m4s]&&UnsameQ[m1s,m3s];
s13s24=SameQ[m1s,m3s]&&SameQ[m2s,m4s]&&UnsameQ[m1s,m2s];
s14s23=SameQ[m1s,m4s]&&SameQ[m2s,m3s]&&UnsameQ[m1s,m2s];
s123=SameQ[m1s,m2s,m3s]&&UnsameQ[m1s,m4s];
s124=SameQ[m1s,m2s,m4s]&&UnsameQ[m1s,m3s];
s134=SameQ[m1s,m3s,m4s]&&UnsameQ[m1s,m2s];
s234=SameQ[m2s,m3s,m4s]&&UnsameQ[m2s,m1s];
s1234=SameQ[m1s,m2s,m3s,m4s];
us=UnsameQ[m1s,m2s,m3s,m4s];
Which[s12,-(t1 (t1+t1^2-2 t3) (-1+t3)^2 Log[t1]-(-1+t1) (-(-1+t3) (-t1^2 (-2+t3)+t1 (-3+t3) t3+t3^2)+(-1+t1)^2 t3^2 Log[t3]))/(2 m4s^2 (-1+t1)^3 (t1-t3)^2 (-1+t3)^2),
s13,-(t1 (t1+t1^2-2 t2) (-1+t2)^2 Log[t1]-(-1+t1) (-(-1+t2) (-t1^2 (-2+t2)+t1 (-3+t2) t2+t2^2)+(-1+t1)^2 t2^2 Log[t2]))/(2 m4s^2 (-1+t1)^3 (t1-t2)^2 (-1+t2)^2),
s14,(t1^2 (-1+t2)^2 (t2-t3) (-1+t3)^2 Log[t1]-(-1+t1) ((-1+t1) t2^2 (t1-t3) (-1+t3)^2 Log[t2]-(t1-t2) (-1+t2) ((t1-t3) (-1+t3) (-t2+t3)+(-1+t1) (-1+t2) t3^2 Log[t3])))/(2 m4s^2 (-1+t1)^2 (t1-t2) (-1+t2)^2 (t1-t3) (t2-t3) (-1+t3)^2),
s23,(t1^2 (-1+t2)^3 Log[t1]-(-1+t1) ((-1+t2) (2 t2^2+t1^2 (1+t2)-t1 t2 (3+t2))+t2 (-2 t1^2-t2 (1+t2)+t1 (2+t2+t2^2)) Log[t2]))/(2 m4s^2 (-1+t1)^2 (t1-t2)^2 (-1+t2)^3),
s24,(2 t1^2 (-1+t3)^3 Log[t1]-(-1+t1) ((-1+t3) (t1+t1^2 (1-3 t3)+3 t1 t3^2-t3 (1+t3))+2 (-1+t1)^2 t3^2 Log[t3]))/(4 m4s^2 (-1+t1)^3 (t1-t3) (-1+t3)^3),
s34,(t1^2 (-1+t2)^3 Log[t1]+(-1+t1) (1/2 (t1-t2) (-1+t2) (-1-t2+t1 (-1+3 t2))-(-1+t1)^2 t2^2 Log[t2]))/(2 m4s^2 (-1+t1)^3 (t1-t2) (-1+t2)^3),
s12s34,-((1+4 t1-5 t1^2+2 t1 (2+t1) Log[t1])/(4 m4s^2 (-1+t1)^4)),
s13s24,-((1+4 t1-5 t1^2+2 t1 (2+t1) Log[t1])/(4 m4s^2 (-1+t1)^4)),
s14s23,-((1+4 t2-5 t2^2+2 t2 (2+t2) Log[t2])/(4 m4s^2 (-1+t2)^4)),
s123,(5-4 t1-t1^2+(2+4 t1) Log[t1])/(4 m4s^2 (-1+t1)^4),
s124,(-1+6 t3-3 t3^2-2 t3^3+6 t3^2 Log[t3])/(12 m4s^2 (-1+t3)^4),
s234,(-1+6 t1-3 t1^2-2 t1^3+6 t1^2 Log[t1])/(12 m4s^2 (-1+t1)^4),
s1234,-(1/(24 m4s^2)),
us,
(t1^2 (-1+t2)^2 (t2-t3) (-1+t3)^2 Log[t1]-(-1+t1) ((-1+t1) t2^2 (t1-t3) (-1+t3)^2 Log[t2]-(t1-t2) (-1+t2) ((t1-t3) (-1+t3) (-t2+t3)+(-1+t1) (-1+t2) t3^2 Log[t3])))/(2 m4s^2 (-1+t1)^2 (t1-t2) (-1+t2)^2 (t1-t3) (t2-t3) (-1+t3)^2)]];


DFant[LoopTools`dd00,m1s_,m2s_,m3s_,m4s_]:=Block[{t1,t2,t3,s12,s13,s14,s23,s24,s34,s12s34,s13s24,s14s23,s123,s124,s134,s234,s1234,us},
{t1,t2,t3}={m1s/m4s,m2s/m4s,m3s/m4s};
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s,m4s];
s13=SameQ[m1s,m3s]&&UnsameQ[m1s,m2s,m4s];
s14=SameQ[m1s,m4s]&&UnsameQ[m1s,m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m2s,m1s,m4s];
s24=SameQ[m2s,m4s]&&UnsameQ[m2s,m1s,m3s];
s34=SameQ[m3s,m4s]&&UnsameQ[m3s,m1s,m2s];
s12s34=SameQ[m1s,m2s]&&SameQ[m3s,m4s]&&UnsameQ[m1s,m3s];
s13s24=SameQ[m1s,m3s]&&SameQ[m2s,m4s]&&UnsameQ[m1s,m2s];
s14s23=SameQ[m1s,m4s]&&SameQ[m2s,m3s]&&UnsameQ[m1s,m2s];
s123=SameQ[m1s,m2s,m3s]&&UnsameQ[m1s,m4s];
s124=SameQ[m1s,m2s,m4s]&&UnsameQ[m1s,m3s];
s134=SameQ[m1s,m3s,m4s]&&UnsameQ[m1s,m2s];
s234=SameQ[m2s,m3s,m4s]&&UnsameQ[m2s,m1s];
s1234=SameQ[m1s,m2s,m3s,m4s];
us=UnsameQ[m1s,m2s,m3s,m4s];
Which[
s12,-(1/2)(-t1 (-1+t3) (t1-2 t3+t1 t3) Log[t1]+(-1+t1) (t1 (t1-t3) (-1+t3)+(-1+t1) t3^2 Log[t3]))/(2 m4s (-1+t1)^2 (t1-t3)^2 (-1+t3)),
s13,-(1/2)(-t1 (-1+t2) (t1-2 t2+t1 t2) Log[t1]+(-1+t1) (t1 (t1-t2) (-1+t2)+(-1+t1) t2^2 Log[t2]))/(2 m4s (-1+t1)^2 (t1-t2)^2 (-1+t2)),
s14,(-t2^2 (-1+t3)^2 Log[t2]+(-1+t2) (-(t2-t3) (-1+t3)+(-1+t2) t3^2 Log[t3]))/(4 m4s (-1+t2)^2 (t2-t3) (-1+t3)^2),
s23,(-t1^2 (-1+t2)^2 Log[t1]+(-1+t1) t2 ((t1-t2) (-1+t2)+(t1 (-2+t2)+t2) Log[t2]))/(4 m4s (-1+t1) (t1-t2)^2 (-1+t2)^2),
s24,(-t1^2 (-1+t3)^2 Log[t1]+(-1+t1) (-(t1-t3) (-1+t3)+(-1+t1) t3^2 Log[t3]))/(4 m4s (-1+t1)^2 (t1-t3) (-1+t3)^2),
s34,(-t1^2 (-1+t2)^2 Log[t1]+(-1+t1) (-(t1-t2) (-1+t2)+(-1+t1) t2^2 Log[t2]))/(4 m4s (-1+t1)^2 (t1-t2) (-1+t2)^2),
s12s34,(1-t1^2+2 t1 Log[t1])/(4 m4s (-1+t1)^3),
s13s24,(1-t1^2+2 t1 Log[t1])/(4 m4s (-1+t1)^3),
s14s23,(1-t2^2+2 t2 Log[t2])/(4 m4s (-1+t2)^3),
s123,-((3-4 t1+t1^2+2 Log[t1])/(8 m4s (-1+t1)^3)),
s124,(1-4 t3+3 t3^2-2 t3^2 Log[t3])/(8 m4s (-1+t3)^3),
s134,(1-4 t2+3 t2^2-2 t2^2 Log[t2])/(8 m4s (-1+t2)^3),
s234,(1-4 t1+3 t1^2-2 t1^2 Log[t1])/(8 m4s (-1+t1)^3),
s1234,-(1/(12 m4s)),
us,-(1/2)(t1^2 (-1+t2) (t2-t3) (-1+t3) Log[t1]-(-1+t1) (t2^2 (t1-t3) (-1+t3) Log[t2]-(t1-t2) (-1+t2) t3^2 Log[t3]))/(2 m4s (-1+t1) (t1-t2) (-1+t2) (t1-t3) (t2-t3) (-1+t3))]
];


DFant[LoopTools`dd11,m1s_,m2s_,m3s_,m4s_]:=Block[{t1,t2,t3,s12,s13,s14,s23,s24,s34,s12s34,s13s24,s14s23,s123,s124,s134,s234,s1234,us},
{t1,t2,t3}={m1s/m4s,m2s/m4s,m3s/m4s};
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s,m4s];
s13=SameQ[m1s,m3s]&&UnsameQ[m1s,m2s,m4s];
s14=SameQ[m1s,m4s]&&UnsameQ[m1s,m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m2s,m1s,m4s];
s24=SameQ[m2s,m4s]&&UnsameQ[m2s,m1s,m3s];
s34=SameQ[m3s,m4s]&&UnsameQ[m3s,m1s,m2s];
s12s34=SameQ[m1s,m2s]&&SameQ[m3s,m4s]&&UnsameQ[m1s,m3s];
s13s24=SameQ[m1s,m3s]&&SameQ[m2s,m4s]&&UnsameQ[m1s,m2s];
s14s23=SameQ[m1s,m4s]&&SameQ[m2s,m3s]&&UnsameQ[m1s,m2s];
s123=SameQ[m1s,m2s,m3s]&&UnsameQ[m1s,m4s];
s124=SameQ[m1s,m2s,m4s]&&UnsameQ[m1s,m3s];
s134=SameQ[m1s,m3s,m4s]&&UnsameQ[m1s,m2s];
s234=SameQ[m2s,m3s,m4s]&&UnsameQ[m2s,m1s];
s1234=SameQ[m1s,m2s,m3s,m4s];
us=UnsameQ[m1s,m2s,m3s,m4s];
Which[
s12,(6 (-1+t3) (6 t1^2 t3^2-t3^3-4 t1^3 t3 (1+t3)+t1^4 (1+t3+t3^2)) Log[t1]+(-1+t1) (t1^5 (-1+t3)+11 (-1+t3) t3^3+t1 t3^2 (18-11 t3-7 t3^2)+t1^4 (5+t3-6 t3^2)+t1^2 t3 (-9+12 t3-5 t3^2+2 t3^3)+t1^3 (2-17 t3+12 t3^2+3 t3^3)-6 (-1+t1)^3 t3^3 Log[t3]))/(18 m4s^2 (-1+t1)^4 (t1-t3)^4 (-1+t3)),
s13,-(-2 t1^2 (-1+t2)^3 (t1^2-3 t2+2 t1 t2) Log[t1]+(-1+t1) ((-1+t2) (t2^3 (1+t2)+t1^2 t2 (3+7 t2-4 t2^2)+t1 t2^2 (-6+t2-t2^2)+t1^3 (2-9 t2+5 t2^2))+2 t2 (-t1 (-4+t2) t2^2-t2^3+t1^3 (3-3 t2+t2^2)+t1^2 (-3+3 t2-5 t2^2+2 t2^3)) Log[t2]))/(6 m4s^2 (-1+t1)^2 (t1-t2)^4 (-1+t2)^3),
s14,(-2 t2 (-1+t3)^2 (t2^2 (1-4 t3)-3 t2 t3+3 t3^2+t2^3 (2+t3)) Log[t2]+(-1+t2) (-(-1+t3) (t2^3 (5-3 t3)-t2^4 (-1+t3)+t2 (11-5 t3) t3^2-2 t3^3+t2^2 t3 (-14+7 t3+t3^2))+2 (-1+t2)^3 t3^3 Log[t3]))/(6 m4s^2 (-1+t2)^4 (t2-t3)^3 (-1+t3)^2),
s23,(-6 t1^3 (-1+t2)^4 Log[t1]+(-1+t1) ((-1+t2) (3 t1 t2^2 (3+5 t2-2 t2^2)+t2^3 (-2-5 t2+t2^2)+3 t1^2 t2 (-6-t2+t2^2)+t1^3 (11-7 t2+2 t2^2))-6 (t1^3-t1 (-4+t2) t2^3-t2^4-t1^2 t2^2 (6-4 t2+t2^2)) Log[t2]))/(18 m4s^2 (-1+t1) (t1-t2)^4 (-1+t2)^4),
s24,(-6 t1^3 (-1+t3)^4 Log[t1]+(-1+t1) ((-1+t3) (t1^3 (-2+7 t3-11 t3^2)+t3 (-1+5 t3+2 t3^2)+t1 (1-12 t3^2-7 t3^3)+t1^2 (-5+12 t3+11 t3^3))+6 (-1+t1)^3 t3^3 Log[t3]))/(18 m4s^2 (-1+t1)^4 (t1-t3) (-1+t3)^4),
s34,(-2 t1^3 (-1+t2)^4 Log[t1]+(-1+t1) ((-1+t2) (t2^3 (5+t2)+t1^2 t2 (11+7 t2)+t1^3 (-2-5 t2+t2^2)-t1 t2^2 (14+3 t2+t2^2))+2 t2 (3 t1^3-t2^2 (1+2 t2)+t1 t2 (3+5 t2+t2^2)+t1^2 (-3-3 t2-4 t2^2+t2^3)) Log[t2]))/(6 m4s^2 (-1+t1)^2 (t1-t2)^3 (-1+t2)^4),
s12s34,(17-9 t1-9 t1^2+t1^3+6 (1+3 t1) Log[t1])/(18 m4s^2 (-1+t1)^5),
s13s24,(-1+9 t1+9 t1^2-17 t1^3+6 t1^2 (3+t1) Log[t1])/(18 m4s^2 (-1+t1)^5),
s14s23,(17-9 t2-9 t2^2+t2^3+6 (1+3 t2) Log[t2])/(18 m4s^2 (-1+t2)^5),
s123,-((5+3/(2 t1)-9 t1+3 t1^2-t1^3/2+6 Log[t1])/(18 m4s^2 (-1+t1)^5)),
s124,(-1+6 t3-18 t3^2+10 t3^3+3 t3^4-12 t3^3 Log[t3])/(36 m4s^2 (-1+t3)^5),
s134,(-1-9 t2+9 t2^2+t2^3-6 t2 (1+t2) Log[t2])/(6 m4s^2 (-1+t2)^5),
s234,(-1+6 t1-18 t1^2+10 t1^3+3 t1^4-12 t1^3 Log[t1])/(36 m4s^2 (-1+t1)^5),
s1234,1/(60 m4s^2),
us,
(-2 t1^3 (-1+t2)^3 (t2-t3)^3 (-1+t3) Log[t1]+(-1+t1) (2 t2 (-1+t3) (t1 t2 (t2^4-6 t2^2 t3^2+3 t2^3 t3^2+3 t3^3-t2 t3^3)-t2^2 t3 (-3 t2^2 t3+t3^2+t2^3 (1+t3))+t1^2 (t2^5+6 t2^3 t3-3 t3^3+3 t2 t3^3-t2^2 t3^3-3 t2^4 (1+t3))+t1^3 (3 t3^2-3 t2 t3 (1+t3)+t2^2 (1+t3+t3^2))) Log[t2]+(t1-t2) (-1+t2) (t2 (t2-t3) (-1+t3) (-t2 t3 (t2+t2^2-3 t3+t2 t3)+t1 (t2^2+t2^3-5 t3^2+3 t2 t3^2)+t1^2 (t2^2+5 t3-3 t2 (1+t3)))-2 (t1-t2)^2 (-1+t2)^2 t3^3 Log[t3])))/(6 m4s^2 (-1+t1) (t1-t2)^3 (-1+t2)^3 (t1-t3) (t2-t3)^3 (-1+t3))]];


DFant[LoopTools`dd12,m1s_,m2s_,m3s_,m4s_]:=Block[{t1,t2,t3,s12,s13,s14,s23,s24,s34,s12s34,s13s24,s14s23,s123,s124,s134,s234,s1234,us},
{t1,t2,t3}={m1s/m4s,m2s/m4s,m3s/m4s};
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s,m4s];
s13=SameQ[m1s,m3s]&&UnsameQ[m1s,m2s,m4s];
s14=SameQ[m1s,m4s]&&UnsameQ[m1s,m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m2s,m1s,m4s];
s24=SameQ[m2s,m4s]&&UnsameQ[m2s,m1s,m3s];
s34=SameQ[m3s,m4s]&&UnsameQ[m3s,m1s,m2s];
s12s34=SameQ[m1s,m2s]&&SameQ[m3s,m4s]&&UnsameQ[m1s,m3s];
s13s24=SameQ[m1s,m3s]&&SameQ[m2s,m4s]&&UnsameQ[m1s,m2s];
s14s23=SameQ[m1s,m4s]&&SameQ[m2s,m3s]&&UnsameQ[m1s,m2s];
s123=SameQ[m1s,m2s,m3s]&&UnsameQ[m1s,m4s];
s124=SameQ[m1s,m2s,m4s]&&UnsameQ[m1s,m3s];
s134=SameQ[m1s,m3s,m4s]&&UnsameQ[m1s,m2s];
s234=SameQ[m2s,m3s,m4s]&&UnsameQ[m2s,m1s];
s1234=SameQ[m1s,m2s,m3s,m4s];
us=UnsameQ[m1s,m2s,m3s,m4s];
Which[
s12,(-2 t1 (-1+t3)^2 (t1^2 (-4+t3) t3+3 t3^2-3 t1 t3^2+t1^3 (1+2 t3)) Log[t1]+(-1+t1) (-(-1+t3) (-t1^4 (-1+t3)+3 t1 (1-3 t3) t3^2+2 t3^3+t1^3 (1+t3-4 t3^2)+t1^2 t3 (-6+7 t3+5 t3^2))+2 (-1+t1)^2 t3^2 (t3^2+t1 (-3+2 t3)) Log[t3]))/(12 m4s^2 (-1+t1)^3 (t1-t3)^4 (-1+t3)^2),
s13,(-2 t1 (-1+t2)^2 (t1^2 (-4+t2) t2+3 t2^2-3 t1 t2^2+t1^3 (1+2 t2)) Log[t1]+(-1+t1) (-(-1+t2) (-t1^4 (-1+t2)+3 t1 (1-3 t2) t2^2+2 t2^3+t1^3 (1+t2-4 t2^2)+t1^2 t2 (-6+7 t2+5 t2^2))+2 (-1+t1)^2 t2^2 (t2^2+t1 (-3+2 t2)) Log[t2]))/(12 m4s^2 (-1+t1)^3 (t1-t2)^4 (-1+t2)^2),
s14,(t2^2 (-1+t3)^3 (t2+t2^2-3 t3+t2 t3) Log[t2]-(-1+t2) (2 (-1+t3) (-t3^3+t2 t3^2 (2+t3)-t2^2 t3 (2+t3^2)+t2^3 (1-t3+t3^2))+(-1+t2)^2 t3^2 (t2 (-3+t3)+t3 (1+t3)) Log[t3]))/(6 m4s^2 (-1+t2)^3 (t2-t3)^3 (-1+t3)^3),
s23,(-6 t1^3 (-1+t2)^4 Log[t1]+(-1+t1) ((-1+t2) (3 t1 t2^2 (3+5 t2-2 t2^2)+t2^3 (-2-5 t2+t2^2)+3 t1^2 t2 (-6-t2+t2^2)+t1^3 (11-7 t2+2 t2^2))-6 (t1^3-t1 (-4+t2) t2^3-t2^4-t1^2 t2^2 (6-4 t2+t2^2)) Log[t2]))/(36 m4s^2 (-1+t1) (t1-t2)^4 (-1+t2)^4),
s24,(-2 t1^3 (-1+t3)^4 Log[t1]+(-1+t1) (-(-1+t3) (t3^2 (1+5 t3)+t1^3 (1-5 t3-2 t3^2)-t1 t3 (2+7 t3+9 t3^2)+t1^2 (1+t3+14 t3^2+2 t3^3))-2 (-1+t1)^2 t3^2 (3 t1-t3 (2+t3)) Log[t3]))/(12 m4s^2 (-1+t1)^3 (t1-t3)^2 (-1+t3)^4),
s34,(-t1^3 (-1+t2)^5 Log[t1]+(-1+t1) (1/2 (t1-t2) (-1+t2)^2 (t2 (1+5 t2)+t1^2 (-1+5 t2+2 t2^2)-t1 (1+2 t2+9 t2^2))+(-1+t1)^2 t2^2 (-3 t1 (-1+t2)+t2 (-2+t2+t2^2)) Log[t2]))/(6 m4s^2 (-1+t1)^3 (t1-t2)^2 (-1+t2)^5),
s12s34,(-1-9 t1+9 t1^2+t1^3-6 t1 (1+t1) Log[t1])/(12 m4s^2 (-1+t1)^5),
s13s24,(-1-9 t1+9 t1^2+t1^3-6 t1 (1+t1) Log[t1])/(12 m4s^2 (-1+t1)^5),
s14s23,(17-9 t2-9 t2^2+t2^3+6 (1+3 t2) Log[t2])/(36 m4s^2 (-1+t2)^5),
s123,-((5+3/(2 t1)-9 t1+3 t1^2-t1^3/2+6 Log[t1])/(36 m4s^2 (-1+t1)^5)),
s124,(-1+9 t3+9 t3^2-17 t3^3+6 t3^2 (3+t3) Log[t3])/(36 m4s^2 (-1+t3)^5),
s134,(-1+9 t2+9 t2^2-17 t2^3+6 t2^2 (3+t2) Log[t2])/(36 m4s^2 (-1+t2)^5),
s234,(-1+6 t1-18 t1^2+10 t1^3+3 t1^4-12 t1^3 Log[t1])/(72 m4s^2 (-1+t1)^5),
s1234,1/(120 m4s^2),
us,(-t1^3 (-1+t2)^2 (t2-t3)^3 (-1+t3)^2 Log[t1]+(-1+t1) (t2^2 (t1-t3)^2 (-1+t3)^2 (t2 (t2^2-2 t3+t2 t3)-t1 (t2-3 t3+2 t2 t3)) Log[t2]+(t1-t2) (-1+t2) ((t1-t3) (t2-t3) (-1+t3) (t2 t3 (t2+t3-2 t2 t3)+t1 (t2^2 (-1+t3)-t3^2+t2 t3^2))+(-1+t2) t3^2 (t2 t3 (t2 (-2+t3)+t3^2)+t1^2 (t3+t2 (-3+2 t3))-t1 (t2 (-1+t3) t3+t3^3+t2^2 (-3+2 t3))) Log[t3])))/(6 m4s^2 (-1+t1) (t1-t2)^2 (-1+t2)^2 (t1-t3)^2 (t2-t3)^3 (-1+t3)^2)]];


DFant[LoopTools`dd13,m1s_,m2s_,m3s_,m4s_]:=Block[{t1,t2,t3,s12,s13,s14,s23,s24,s34,s12s34,s13s24,s14s23,s123,s124,s134,s234,s1234,us},
{t1,t2,t3}={m1s/m4s,m2s/m4s,m3s/m4s};
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s,m4s];
s13=SameQ[m1s,m3s]&&UnsameQ[m1s,m2s,m4s];
s14=SameQ[m1s,m4s]&&UnsameQ[m1s,m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m2s,m1s,m4s];
s24=SameQ[m2s,m4s]&&UnsameQ[m2s,m1s,m3s];
s34=SameQ[m3s,m4s]&&UnsameQ[m3s,m1s,m2s];
s12s34=SameQ[m1s,m2s]&&SameQ[m3s,m4s]&&UnsameQ[m1s,m3s];
s13s24=SameQ[m1s,m3s]&&SameQ[m2s,m4s]&&UnsameQ[m1s,m2s];
s14s23=SameQ[m1s,m4s]&&SameQ[m2s,m3s]&&UnsameQ[m1s,m2s];
s123=SameQ[m1s,m2s,m3s]&&UnsameQ[m1s,m4s];
s124=SameQ[m1s,m2s,m4s]&&UnsameQ[m1s,m3s];
s134=SameQ[m1s,m3s,m4s]&&UnsameQ[m1s,m2s];
s234=SameQ[m2s,m3s,m4s]&&UnsameQ[m2s,m1s];
s1234=SameQ[m1s,m2s,m3s,m4s];
us=UnsameQ[m1s,m2s,m3s,m4s];
Which[
s12,(-2 t1 (-1+t3)^2 (t1^2 (1-4 t3)-3 t1 t3+3 t3^2+t1^3 (2+t3)) Log[t1]+(-1+t1) (-(-1+t3) (t1^3 (5-3 t3)-t1^4 (-1+t3)+t1 (11-5 t3) t3^2-2 t3^3+t1^2 t3 (-14+7 t3+t3^2))+2 (-1+t1)^3 t3^3 Log[t3]))/(12 m4s^2 (-1+t1)^4 (t1-t3)^3 (-1+t3)^2),
s13,-(-t1^2 (-1+t2)^3 (t1+t1^2-3 t2+t1 t2) Log[t1]+(-1+t1) (2 (-1+t2) (-t2^3+t1 t2^2 (2+t2)-t1^2 t2 (2+t2^2)+t1^3 (1-t2+t2^2))+(-1+t1)^2 t2^2 (t1 (-3+t2)+t2 (1+t2)) Log[t2]))/(6 m4s^2 (-1+t1)^3 (t1-t2)^3 (-1+t2)^3),
s14,(2 t2^2 (2 t2+t2^2-3 t3) (-1+t3)^3 Log[t2]-(-1+t2) ((-1+t3) (t3^2 (1+t3)+t2 t3 (-2+t3-5 t3^2)+t2^3 (5-9 t3+2 t3^2)+t2^2 (1-7 t3+14 t3^2-2 t3^3))+2 (-1+t2)^3 t3^3 Log[t3]))/(12 m4s^2 (-1+t2)^4 (t2-t3)^2 (-1+t3)^3),
s23,(-2 t1^3 (-1+t2)^4 Log[t1]+(-1+t1) ((-1+t2) (t2^3 (5+t2)+t1^2 t2 (11+7 t2)+t1^3 (-2-5 t2+t2^2)-t1 t2^2 (14+3 t2+t2^2))+2 t2 (3 t1^3-t2^2 (1+2 t2)+t1 t2 (3+5 t2+t2^2)+t1^2 (-3-3 t2-4 t2^2+t2^3)) Log[t2]))/(12 m4s^2 (-1+t1)^2 (t1-t2)^3 (-1+t2)^4),
s24,(-6 t1^3 (-1+t3)^4 Log[t1]+(-1+t1) ((-1+t3) (t1^3 (-2+7 t3-11 t3^2)+t3 (-1+5 t3+2 t3^2)+t1 (1-12 t3^2-7 t3^3)+t1^2 (-5+12 t3+11 t3^3))+6 (-1+t1)^3 t3^3 Log[t3]))/(36 m4s^2 (-1+t1)^4 (t1-t3) (-1+t3)^4),
s34,(-2 t1^3 (-1+t2)^4 Log[t1]+(-1+t1) (-(-1+t2) (t2^2 (1+5 t2)+t1^3 (1-5 t2-2 t2^2)-t1 t2 (2+7 t2+9 t2^2)+t1^2 (1+t2+14 t2^2+2 t2^3))-2 (-1+t1)^2 t2^2 (3 t1-t2 (2+t2)) Log[t2]))/(12 m4s^2 (-1+t1)^3 (t1-t2)^2 (-1+t2)^4),
s12s34,(-1-9 t1+9 t1^2+t1^3-6 t1 (1+t1) Log[t1])/(12 m4s^2 (-1+t1)^5),
s13s24,(-1+9 t1+9 t1^2-17 t1^3+6 t1^2 (3+t1) Log[t1])/(36 m4s^2 (-1+t1)^5),
s14s23,(-1-9 t2+9 t2^2+t2^3-6 t2 (1+t2) Log[t2])/(12 m4s^2 (-1+t2)^5),
s123,(17-9 t1-9 t1^2+t1^3+6 (1+3 t1) Log[t1])/(36 m4s^2 (-1+t1)^5),
s124,(-1+6 t3-18 t3^2+10 t3^3+3 t3^4-12 t3^3 Log[t3])/(72 m4s^2 (-1+t3)^5),
s134,(-1+9 t2+9 t2^2-17 t2^3+6 t2^2 (3+t2) Log[t2])/(36 m4s^2 (-1+t2)^5),
s234,(-1+6 t1-18 t1^2+10 t1^3+3 t1^4-12 t1^3 Log[t1])/(72 m4s^2 (-1+t1)^5),
s1234,1/(120 m4s^2),
us,
(-t1^3 (-1+t2)^3 (t2-t3)^2 (-1+t3)^2 Log[t1]+(-1+t1) (-(-1+t1) t2^2 (-1+t3)^2 (t2 (t2+t2^2-2 t3) t3-t1 (t2^2+t2^3-3 t3^2+t2 t3^2)+t1^2 (-3 t3+t2 (2+t3))) Log[t2]+(t1-t2) (-1+t2) ((t1-t3) (t2-t3) (-1+t3) (-t2 (t2 (-2+t3)+t3)+t1 (-t2+t2^2 (-1+t3)+t3))+(-1+t1) (t1-t2) (-1+t2)^2 t3^3 Log[t3])))/(6 m4s^2 (-1+t1)^2 (t1-t2)^2 (-1+t2)^3 (t1-t3) (t2-t3)^2 (-1+t3)^2)]];


DFant[LoopTools`dd22,m1s_,m2s_,m3s_,m4s_]:=Block[{t1,t2,t3,s12,s13,s14,s23,s24,s34,s12s34,s13s24,s14s23,s123,s124,s134,s234,s1234,us},
{t1,t2,t3}={m1s/m4s,m2s/m4s,m3s/m4s};
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s,m4s];
s13=SameQ[m1s,m3s]&&UnsameQ[m1s,m2s,m4s];
s14=SameQ[m1s,m4s]&&UnsameQ[m1s,m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m2s,m1s,m4s];
s24=SameQ[m2s,m4s]&&UnsameQ[m2s,m1s,m3s];
s34=SameQ[m3s,m4s]&&UnsameQ[m3s,m1s,m2s];
s12s34=SameQ[m1s,m2s]&&SameQ[m3s,m4s]&&UnsameQ[m1s,m3s];
s13s24=SameQ[m1s,m3s]&&SameQ[m2s,m4s]&&UnsameQ[m1s,m2s];
s14s23=SameQ[m1s,m4s]&&SameQ[m2s,m3s]&&UnsameQ[m1s,m2s];
s123=SameQ[m1s,m2s,m3s]&&UnsameQ[m1s,m4s];
s124=SameQ[m1s,m2s,m4s]&&UnsameQ[m1s,m3s];
s134=SameQ[m1s,m3s,m4s]&&UnsameQ[m1s,m2s];
s234=SameQ[m2s,m3s,m4s]&&UnsameQ[m2s,m1s];
s1234=SameQ[m1s,m2s,m3s,m4s];
us=UnsameQ[m1s,m2s,m3s,m4s];
Which[
s12,-(-2 t1^2 (-1+t3)^3 (t1^2-3 t3+2 t1 t3) Log[t1]+(-1+t1) ((-1+t3) (t3^3 (1+t3)+t1^2 t3 (3+7 t3-4 t3^2)+t1 t3^2 (-6+t3-t3^2)+t1^3 (2-9 t3+5 t3^2))+2 t3 (-t1 (-4+t3) t3^2-t3^3+t1^3 (3-3 t3+t3^2)+t1^2 (-3+3 t3-5 t3^2+2 t3^3)) Log[t3]))/(6 m4s^2 (-1+t1)^2 (t1-t3)^4 (-1+t3)^3),
s13,(6 (-1+t2) (6 t1^2 t2^2-t2^3-4 t1^3 t2 (1+t2)+t1^4 (1+t2+t2^2)) Log[t1]+(-1+t1) (t1^5 (-1+t2)+11 (-1+t2) t2^3+t1 t2^2 (18-11 t2-7 t2^2)+t1^4 (5+t2-6 t2^2)+t1^2 t2 (-9+12 t2-5 t2^2+2 t2^3)+t1^3 (2-17 t2+12 t2^2+3 t2^3)-6 (-1+t1)^3 t2^3 Log[t2]))/(18 m4s^2 (-1+t1)^4 (t1-t2)^4 (-1+t2)),
s14,(-2 t2^3 (-1+t3)^4 Log[t2]+(-1+t2) ((-1+t3) (t3^3 (5+t3)+t2^2 t3 (11+7 t3)+t2^3 (-2-5 t3+t3^2)-t2 t3^2 (14+3 t3+t3^2))+2 t3 (3 t2^3-t3^2 (1+2 t3)+t2 t3 (3+5 t3+t3^2)+t2^2 (-3-3 t3-4 t3^2+t3^3)) Log[t3]))/(6 m4s^2 (-1+t2)^2 (t2-t3)^3 (-1+t3)^4),
s23,(-6 t1^3 (-1+t2)^4 Log[t1]+(-1+t1) ((-1+t2) (3 t1 t2^2 (3+5 t2-2 t2^2)+t2^3 (-2-5 t2+t2^2)+3 t1^2 t2 (-6-t2+t2^2)+t1^3 (11-7 t2+2 t2^2))-6 (t1^3-t1 (-4+t2) t2^3-t2^4-t1^2 t2^2 (6-4 t2+t2^2)) Log[t2]))/(18 m4s^2 (-1+t1) (t1-t2)^4 (-1+t2)^4),
s24,(-2 t1^3 (-1+t3)^4 Log[t1]+(-1+t1) ((-1+t3) (t3^3 (5+t3)+t1^2 t3 (11+7 t3)+t1^3 (-2-5 t3+t3^2)-t1 t3^2 (14+3 t3+t3^2))+2 t3 (3 t1^3-t3^2 (1+2 t3)+t1 t3 (3+5 t3+t3^2)+t1^2 (-3-3 t3-4 t3^2+t3^3)) Log[t3]))/(6 m4s^2 (-1+t1)^2 (t1-t3)^3 (-1+t3)^4),
s34,(-2 t1^3 (-1+t2)^4 Log[t1]+(-1+t1) (-(1/3) (t1-t2) (-1+t2) (-1+5 t2+2 t2^2+t1 (5-10 t2-7 t2^2)+t1^2 (2-7 t2+11 t2^2))+2 (-1+t1)^3 t2^3 Log[t2]))/(6 m4s^2 (-1+t1)^4 (t1-t2) (-1+t2)^4),
s12s34,(-1+9 t1+9 t1^2-17 t1^3+6 t1^2 (3+t1) Log[t1])/(18 m4s^2 (-1+t1)^5),
s13s24,(17-9 t1-9 t1^2+t1^3+6 (1+3 t1) Log[t1])/(18 m4s^2 (-1+t1)^5),
s14s23,(17-9 t2-9 t2^2+t2^3+6 (1+3 t2) Log[t2])/(18 m4s^2 (-1+t2)^5),
s123,-((5+3/(2 t1)-9 t1+3 t1^2-t1^3/2+6 Log[t1])/(18 m4s^2 (-1+t1)^5)),
s124,(-1-9 t3+9 t3^2+t3^3-6 t3 (1+t3) Log[t3])/(6 m4s^2 (-1+t3)^5),
s134,(-1+6 t2-18 t2^2+10 t2^3+3 t2^4-12 t2^3 Log[t2])/(36 m4s^2 (-1+t2)^5),
s234,(-1+6 t1-18 t1^2+10 t1^3+3 t1^4-12 t1^3 Log[t1])/(36 m4s^2 (-1+t1)^5),
s1234,1/(60 m4s^2),
us,(-2 t1^3 (-1+t2) (t2-t3)^3 (-1+t3)^3 Log[t1]+(-1+t1) (2 t2^3 (t1-t3)^3 (-1+t3)^3 Log[t2]-(t1-t2) (-1+t2) t3 ((t1-t3) (-1+t3) (-t2+t3) (t1 (t2 (5-3 t3)+(-3+t3) t3)+t3 (t2 (-3+t3)+t3 (1+t3)))+2 (t3^2 (t2^2+t2 (-3+t3) t3^2+t3^3)+t1^2 (t2 (-3+t3) t3+t3^2+t2^2 (3-3 t3+t3^2))+t1 t3 (t2^2 (-3+t3)+(-3+t3) t3^3+t2 (t3+6 t3^2-3 t3^3))) Log[t3])))/(6 m4s^2 (-1+t1) (t1-t2) (-1+t2) (t1-t3)^3 (t2-t3)^3 (-1+t3)^3)]];


DFant[LoopTools`dd23,m1s_,m2s_,m3s_,m4s_]:=Block[{t1,t2,t3,s12,s13,s14,s23,s24,s34,s12s34,s13s24,s14s23,s123,s124,s134,s234,s1234,us},
{t1,t2,t3}={m1s/m4s,m2s/m4s,m3s/m4s};
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s,m4s];
s13=SameQ[m1s,m3s]&&UnsameQ[m1s,m2s,m4s];
s14=SameQ[m1s,m4s]&&UnsameQ[m1s,m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m2s,m1s,m4s];
s24=SameQ[m2s,m4s]&&UnsameQ[m2s,m1s,m3s];
s34=SameQ[m3s,m4s]&&UnsameQ[m3s,m1s,m2s];
s12s34=SameQ[m1s,m2s]&&SameQ[m3s,m4s]&&UnsameQ[m1s,m3s];
s13s24=SameQ[m1s,m3s]&&SameQ[m2s,m4s]&&UnsameQ[m1s,m2s];
s14s23=SameQ[m1s,m4s]&&SameQ[m2s,m3s]&&UnsameQ[m1s,m2s];
s123=SameQ[m1s,m2s,m3s]&&UnsameQ[m1s,m4s];
s124=SameQ[m1s,m2s,m4s]&&UnsameQ[m1s,m3s];
s134=SameQ[m1s,m3s,m4s]&&UnsameQ[m1s,m2s];
s234=SameQ[m2s,m3s,m4s]&&UnsameQ[m2s,m1s];
s1234=SameQ[m1s,m2s,m3s,m4s];
us=UnsameQ[m1s,m2s,m3s,m4s];
Which[
s12,-(-t1^2 (-1+t3)^3 (t1+t1^2-3 t3+t1 t3) Log[t1]+(-1+t1) (2 (-1+t3) (-t3^3+t1 t3^2 (2+t3)-t1^2 t3 (2+t3^2)+t1^3 (1-t3+t3^2))+(-1+t1)^2 t3^2 (t1 (-3+t3)+t3 (1+t3)) Log[t3]))/(6 m4s^2 (-1+t1)^3 (t1-t3)^3 (-1+t3)^3),
s13,(-2 t1 (-1+t2)^2 (t1^2 (1-4 t2)-3 t1 t2+3 t2^2+t1^3 (2+t2)) Log[t1]+(-1+t1) (-(-1+t2) (t1^3 (5-3 t2)-t1^4 (-1+t2)+t1 (11-5 t2) t2^2-2 t2^3+t1^2 t2 (-14+7 t2+t2^2))+2 (-1+t1)^3 t2^3 Log[t2]))/(12 m4s^2 (-1+t1)^4 (t1-t2)^3 (-1+t2)^2),
s14,(-2 t2^3 (-1+t3)^4 Log[t2]+(-1+t2) (-(-1+t3) (t3^2 (1+5 t3)+t2^3 (1-5 t3-2 t3^2)-t2 t3 (2+7 t3+9 t3^2)+t2^2 (1+t3+14 t3^2+2 t3^3))-2 (-1+t2)^2 t3^2 (3 t2-t3 (2+t3)) Log[t3]))/(12 m4s^2 (-1+t2)^3 (t2-t3)^2 (-1+t3)^4),
s23,(-2 t1^3 (-1+t2)^4 Log[t1]+(-1+t1) ((-1+t2) (t2^3 (5+t2)+t1^2 t2 (11+7 t2)+t1^3 (-2-5 t2+t2^2)-t1 t2^2 (14+3 t2+t2^2))+2 t2 (3 t1^3-t2^2 (1+2 t2)+t1 t2 (3+5 t2+t2^2)+t1^2 (-3-3 t2-4 t2^2+t2^3)) Log[t2]))/(12 m4s^2 (-1+t1)^2 (t1-t2)^3 (-1+t2)^4),
s24,(-2 t1^3 (-1+t3)^4 Log[t1]+(-1+t1) (-(-1+t3) (t3^2 (1+5 t3)+t1^3 (1-5 t3-2 t3^2)-t1 t3 (2+7 t3+9 t3^2)+t1^2 (1+t3+14 t3^2+2 t3^3))-2 (-1+t1)^2 t3^2 (3 t1-t3 (2+t3)) Log[t3]))/(12 m4s^2 (-1+t1)^3 (t1-t3)^2 (-1+t3)^4),
s34,(-6 t1^3 (-1+t2)^4 Log[t1]+(-1+t1) ((-1+t2) (t1^3 (-2+7 t2-11 t2^2)+t2 (-1+5 t2+2 t2^2)+t1 (1-12 t2^2-7 t2^3)+t1^2 (-5+12 t2+11 t2^3))+6 (-1+t1)^3 t2^3 Log[t2]))/(36 m4s^2 (-1+t1)^4 (t1-t2) (-1+t2)^4),
s12s34,(-1+9 t1+9 t1^2-17 t1^3+6 t1^2 (3+t1) Log[t1])/(36 m4s^2 (-1+t1)^5),
s13s24,(-1-9 t1+9 t1^2+t1^3-6 t1 (1+t1) Log[t1])/(12 m4s^2 (-1+t1)^5),
s14s23,(-1-9 t2+9 t2^2+t2^3-6 t2 (1+t2) Log[t2])/(12 m4s^2 (-1+t2)^5),
s123,(17-9 t1-9 t1^2+t1^3+6 (1+3 t1) Log[t1])/(36 m4s^2 (-1+t1)^5),
s124,(-1+9 t3+9 t3^2-17 t3^3+6 t3^2 (3+t3) Log[t3])/(36 m4s^2 (-1+t3)^5),
s134,(-1+6 t2-18 t2^2+10 t2^3+3 t2^4-12 t2^3 Log[t2])/(72 m4s^2 (-1+t2)^5),
s234,(-1+6 t1-18 t1^2+10 t1^3+3 t1^4-12 t1^3 Log[t1])/(72 m4s^2 (-1+t1)^5),
s1234,1/(120 m4s^2),
us,
(-t1^3 (-1+t2)^2 (t2-t3)^2 (-1+t3)^3 Log[t1]+(-1+t1) ((-1+t1) t2^3 (t1-t3)^2 (-1+t3)^3 Log[t2]-(t1-t2) (-1+t2) ((t1-t3) (t2-t3) (-1+t3) (-t3 (t2-2 t3+t2 t3)+t1 (-t3 (1+t3)+t2 (1+t3^2)))+(-1+t1) (-1+t2) t3^2 (t1 (t2 (-3+t3)+2 t3)-t3 (-2 t2+t3+t3^2)) Log[t3])))/(6 m4s^2 (-1+t1)^2 (t1-t2) (-1+t2)^2 (t1-t3)^2 (t2-t3)^2 (-1+t3)^3)]];


DFant[LoopTools`dd33,m1s_,m2s_,m3s_,m4s_]:=Block[{t1,t2,t3,s12,s13,s14,s23,s24,s34,s12s34,s13s24,s14s23,s123,s124,s134,s234,s1234,us},
{t1,t2,t3}={m1s/m4s,m2s/m4s,m3s/m4s};
s12=SameQ[m1s,m2s]&&UnsameQ[m1s,m3s,m4s];
s13=SameQ[m1s,m3s]&&UnsameQ[m1s,m2s,m4s];
s14=SameQ[m1s,m4s]&&UnsameQ[m1s,m2s,m3s];
s23=SameQ[m2s,m3s]&&UnsameQ[m2s,m1s,m4s];
s24=SameQ[m2s,m4s]&&UnsameQ[m2s,m1s,m3s];
s34=SameQ[m3s,m4s]&&UnsameQ[m3s,m1s,m2s];
s12s34=SameQ[m1s,m2s]&&SameQ[m3s,m4s]&&UnsameQ[m1s,m3s];
s13s24=SameQ[m1s,m3s]&&SameQ[m2s,m4s]&&UnsameQ[m1s,m2s];
s14s23=SameQ[m1s,m4s]&&SameQ[m2s,m3s]&&UnsameQ[m1s,m2s];
s123=SameQ[m1s,m2s,m3s]&&UnsameQ[m1s,m4s];
s124=SameQ[m1s,m2s,m4s]&&UnsameQ[m1s,m3s];
s134=SameQ[m1s,m3s,m4s]&&UnsameQ[m1s,m2s];
s234=SameQ[m2s,m3s,m4s]&&UnsameQ[m2s,m1s];
s1234=SameQ[m1s,m2s,m3s,m4s];
us=UnsameQ[m1s,m2s,m3s,m4s];
Which[
s12,(2 t1^2 (2 t1+t1^2-3 t3) (-1+t3)^3 Log[t1]-(-1+t1) ((-1+t3) (t3^2 (1+t3)+t1 t3 (-2+t3-5 t3^2)+t1^3 (5-9 t3+2 t3^2)+t1^2 (1-7 t3+14 t3^2-2 t3^3))+2 (-1+t1)^3 t3^3 Log[t3]))/(6 m4s^2 (-1+t1)^4 (t1-t3)^2 (-1+t3)^3),
s13,(2 t1^2 (2 t1+t1^2-3 t2) (-1+t2)^3 Log[t1]-(-1+t1) ((-1+t2) (t2^2 (1+t2)+t1 t2 (-2+t2-5 t2^2)+t1^3 (5-9 t2+2 t2^2)+t1^2 (1-7 t2+14 t2^2-2 t2^3))+2 (-1+t1)^3 t2^3 Log[t2]))/(6 m4s^2 (-1+t1)^4 (t1-t2)^2 (-1+t2)^3),
s14,(-6 t2^3 (-1+t3)^4 Log[t2]+(-1+t2) ((-1+t3) (t2^3 (-2+7 t3-11 t3^2)+t3 (-1+5 t3+2 t3^2)+t2 (1-12 t3^2-7 t3^3)+t2^2 (-5+12 t3+11 t3^3))+6 (-1+t2)^3 t3^3 Log[t3]))/(18 m4s^2 (-1+t2)^4 (t2-t3) (-1+t3)^4),
s23,(-2 t1^3 (-1+t2)^4 Log[t1]+(-1+t1) (-(-1+t2) (t2^2 (1+5 t2)+t1^3 (1-5 t2-2 t2^2)-t1 t2 (2+7 t2+9 t2^2)+t1^2 (1+t2+14 t2^2+2 t2^3))-2 (-1+t1)^2 t2^2 (3 t1-t2 (2+t2)) Log[t2]))/(6 m4s^2 (-1+t1)^3 (t1-t2)^2 (-1+t2)^4),
s24,(-6 t1^3 (-1+t3)^4 Log[t1]+(-1+t1) ((-1+t3) (t1^3 (-2+7 t3-11 t3^2)+t3 (-1+5 t3+2 t3^2)+t1 (1-12 t3^2-7 t3^3)+t1^2 (-5+12 t3+11 t3^3))+6 (-1+t1)^3 t3^3 Log[t3]))/(18 m4s^2 (-1+t1)^4 (t1-t3) (-1+t3)^4),
s34,
(-2 t1^3 (-1+t2)^4 Log[t1]+(-1+t1) (-(1/3) (t1-t2) (-1+t2) (-1+5 t2+2 t2^2+t1 (5-10 t2-7 t2^2)+t1^2 (2-7 t2+11 t2^2))+2 (-1+t1)^3 t2^3 Log[t2]))/(6 m4s^2 (-1+t1)^4 (t1-t2) (-1+t2)^4),
s12s34,(-1+9 t1+9 t1^2-17 t1^3+6 t1^2 (3+t1) Log[t1])/(18 m4s^2 (-1+t1)^5),
s13s24,(-1+9 t1+9 t1^2-17 t1^3+6 t1^2 (3+t1) Log[t1])/(18 m4s^2 (-1+t1)^5),
s14s23,(-1+9 t2+9 t2^2-17 t2^3+6 t2^2 (3+t2) Log[t2])/(18 m4s^2 (-1+t2)^5),
s123,(-1-9 t1+9 t1^2+t1^3-6 t1 (1+t1) Log[t1])/(6 m4s^2 (-1+t1)^5),
s124,(-1+6 t3-18 t3^2+10 t3^3+3 t3^4-12 t3^3 Log[t3])/(36 m4s^2 (-1+t3)^5),
s134,(-1+6 t2-18 t2^2+10 t2^3+3 t2^4-12 t2^3 Log[t2])/(36 m4s^2 (-1+t2)^5),
s234,(-1+6 t1-18 t1^2+10 t1^3+3 t1^4-12 t1^3 Log[t1])/(36 m4s^2 (-1+t1)^5),
s1234,1/(60 m4s^2),
us,
(-2 t1^3 (-1+t2)^3 (t2-t3) (-1+t3)^3 Log[t1]+(-1+t1) (2 (-1+t1)^2 t2^3 (t1-t3) (-1+t3)^3 Log[t2]-(t1-t2) (-1+t2) ((t1-t3) (-1+t3) (-t2+t3) (1+t2+t3-3 t2 t3+t1 (1-3 t3+t2 (-3+5 t3)))+2 (-1+t1)^2 (-1+t2)^2 t3^3 Log[t3])))/(6 m4s^2 (-1+t1)^3 (t1-t2) (-1+t2)^3 (t1-t3) (t2-t3) (-1+t3)^3)]];


DFDant[id_, p1s_, p2s_, p3s_, p4s_, p12s_, p23s_, m1s_, m2s_, m3s_, m4s_] := Block[{}, Print["Derivative for D0i not implemented"]; D0i[id, p1s, p2s, p3s, p4s, p12s, p23s, m1s, m3s, m3s, m4s] - DFant[id, m1s, m2s, m3s, m4s]];


(* ::Subsection:: *)
(*The Main Functions*)


ANT[expr_]:=ReleaseHold[expr/.{A0:>A0ant,B0i:>B0ant,C0i:>C0ant,D0i:>D0ant}];


A0ant[ms_]:=AFant[ms];


B0ant[id_,ps_,m1s_,m2s_]:=Block[{mas,mom,tmas,rmas,mas0,mas1,expr,exprD,x1,x2},
If[(id===LoopTools`bb00||id===LoopTools`bb11)\[And]ps=!=0,
Print["Derivative for B0i[",id,"] not implemented"];
Return[B0i[id, ps, m1s, m2s]];
];
mas={m1s,m2s};
tmas={x1,x2};
rmas=MapThread[{#1,#2}&,{tmas,mas}];
mas0=Map[#[[1]]->#[[2]]&,Cases[rmas,{_,0}]];
mas1=Map[#[[1]]->#[[2]]&,DeleteCases[rmas,{_,0}]];
tmas=tmas/.mas1;
expr=Apply[BFant,Flatten[{id,tmas}]];
expr=Fold[Limit,expr,mas0];
If[ps=!=0,
exprD=Apply[BFDant,Flatten[{id,ps,tmas}]];
expr=expr+Fold[Limit,exprD,mas0];
];
Return[expr];
];


C0ant[id_,p1s_,p2s_,p12s_,m1s_,m2s_,m3s_]:=Block[{mas,mom,tmas,mas0,mas1,expr,exprD,x1,x2,x3,talMas},
If[Not[p1s===0\[And]p2s===0\[And]p12s===0]\[And](id===LoopTools`cc001||id===LoopTools`cc002||id===LoopTools`cc111||id===LoopTools`cc112||id===LoopTools`cc121||id===LoopTools`cc211||id===LoopTools`cc122||id===LoopTools`cc221||id===LoopTools`cc212||id===LoopTools`cc222),
Print["Derivative for C0i[",id,"] not implemented"];
Return[C0i[id, p1s, p2s, p12s, m1s, m2s, m3s]]
];
mas={m1s,m2s,m3s};
tmas={x1,x2,x3};
talMas=Tally[mas];
If[Length[talMas]==3,
Block[{rmas},
rmas=MapThread[{#1,#2}&,{tmas,mas}];
mas0=Map[#[[1]]->#[[2]]&,Cases[rmas,{_,0}]];
mas1=Map[#[[1]]->#[[2]]&,DeleteCases[rmas,{_,0}]]],
If[Length[talMas]==1,
mas0=MapThread[#1->#2&,{tmas,mas}];
mas1={};
];
If[Length[talMas]==2,
Block[{temp,pos,rmas},
temp=Select[talMas, #[[2]] == 2 &][[1,1]];
pos=Flatten[Position[mas,temp]];
rmas=MapThread[{#1,#2}&,{tmas,mas}];
mas0=Map[#[[1]]->#[[2]]&,rmas[[pos]]];
rmas=rmas[[Complement[Range[3],pos]]];
mas0=Join[mas0,Map[#[[1]]->#[[2]]&,Cases[rmas,{_,0}]]];
mas1=Map[#[[1]]->#[[2]]&,DeleteCases[rmas,{_,0}]]];
]];
tmas=tmas/.mas1;
expr=Apply[CFant,Flatten[{id,tmas}]];
expr=Fold[Limit,expr,mas0];
If[Not[p1s===0\[And]p2s===0\[And]p12s===0],
exprD=Apply[CFDant,Flatten[{id,{p1s,p2s,p12s},tmas}]];
expr=expr+Fold[Limit,exprD,mas0];
];
Return[expr];
];
(*C0ant[id_,p1s_,p2s_,p12s_,m1s_,m2s_,m3s_]:=Block[{mas,mom,tmas,rmas,mas0,mas1,expr,exprD,x1,x2,x3},
If[(id===LoopTools`cc22||id===LoopTools`cc12||id===LoopTools`cc11)\[And]Not[p1s===0\[And]p2s===0\[And]p12s===0],
Print["Derivative for C0i[",id,"] not implemented"];
Return[C0i[id, p1s, p2s, p12s, m1s, m2s, m3s]]
];
mas={m1s,m2s,m3s};
tmas={x1,x2,x3};
rmas=MapThread[{#1,#2}&,{tmas,mas}];
mas0=Map[#[[1]]\[Rule]#[[2]]&,Cases[rmas,{_,0}]];
mas1=Map[#[[1]]\[Rule]#[[2]]&,DeleteCases[rmas,{_,0}]];
tmas=tmas/.mas1;
expr=Apply[CFant,Flatten[{id,tmas}]];
expr=Fold[Limit,expr,mas0];
If[Not[p1s===0\[And]p2s===0\[And]p12s===0],
exprD=Apply[CFDant,Flatten[{id,{p1s,p2s,p12s},tmas}]];
expr=expr+Fold[Limit,exprD,mas0];
];
Return[expr];
];*)


D0ant[id_,p1s_,p2s_,p3s_,p4s_,p12s_,p23s_,m1s_,m2s_,m3s_,m4s_]:=Block[{mas,mom,tmas,rmas,mas0,mas1,expr,exprD,x1,x2,x3,x4},
If[Not[p1s===0\[And]p2s===0\[And]p3s===0\[And]p4s===0\[And]p12s===0\[And]p23s===0],
Print["Derivative for D0i not implemented"];
Return[D0i[id, p1s, p2s, p3s, p4s, p12s, p23s, m1s, m3s, m3s, m4s]]
];
mas={m1s,m2s,m3s,m4s};
tmas={x1,x2,x3,x4};
rmas=MapThread[{#1,#2}&,{tmas,mas}];
mas0=Map[#[[1]]->#[[2]]&,Cases[rmas,{_,0}]];
mas1=Map[#[[1]]->#[[2]]&,DeleteCases[rmas,{_,0}]];
tmas=tmas/.mas1;
expr=Apply[DFant,Flatten[{id,tmas}]];
expr=Fold[Limit,expr,mas0];
(*exprD=Apply[`Private`DFDant,Flatten[{id,{p1s,p2s,p3s,p4s,p12s,p23s},tmas}]];
exprD=Fold[Limit,exprD,mas0];*)
Return[expr];
];
