// Command ' <"mol.s" 's, to execute.
// Pi-orbital model on an hexagon.
LIB "primdec.lib";
option(redSB);
proc molsolver()
{
    int n=12;
ring r=0,(y(1..n),e),lp;
//        9      10
//    1      4       7
//       3       6
//    2      5       8
//        11    12
matrix v[12][2]=1,3,2,3,3,4,3,5,4,6,5,6,6,7,6,8,9,4,10,4,11,5,12,5;
matrix p[n][1];
matrix HH[n][n];
int i,i1,i2;
for (i=1; i<=n;i++)
{
p[i,1]=y(i);
}
//print(p)；
print(v);
int k1,k2;
for (i=1;i<=12;i++)
{
 k1=int(v[i,1]);
 k2=int(v[i,2]);
// print(k1);print(k2);
 HH[k1,k2]=-1;
 HH[k2,k1]=-1;
}
for (i=4;i<=6;i++)
{
 k1=int(v[i,1]);
 k2=int(v[i,2]);
// print(k1);print(k2);
 HH[k1,k2]=-2;
 HH[k2,k1]=-2;
}

matrix MI[n][n];
for(i=1;i<=12;i++){ MI[i,i]=1;}
print(HH);
print(det(HH-e*MI));
ideal J=HH*p-e*p,transpose(p)*p-1;
//Krull-dimension is not zero, owing to the degeneracy.
print("Defining Ideal");
print(J);
print("Standard Basis");
std(J);
//The first entry of the standard bases is the polynomial equation for the energy e. What is the galois group?
//Apply primary ideal decomposition. Check the dimension for each of the decomposed components.
list A=primdecGTZ(J);
print("Primary Ideal Decomposition");
print(A);
for(int s=1;s<=size(A);s++)
{
    ideal ps=A[s][2];
    print(std(ps));
    print(dim(std(ps)));
}
return(0);
}

//molsolver(5);
//molsolver(6);


molsolver()
quit;

//For N=5, we obtain a polynomial for y(5). The Galois Group of this polynomial is computed by GAP. 
//gap> y:=Indeterminate(Rationals,1);
//y
//gap> GaloisType(4*y^4+4*y^3-6*y^2-6*y-1);
//1
//gap> TransitiveGroup(4,1);
//C(4) = 4
